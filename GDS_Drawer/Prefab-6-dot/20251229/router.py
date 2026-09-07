import gdstk
import numpy as np
import math

class Router:
    """
    一个全新的布线管理器，负责全局规划和DRC检查。
    """
    def __init__(self, chip_center, active_region_bbox, quantum_dot_polygons, layer, pad_connection_extension_len=10.0):
        self.chip_center = chip_center
        self.active_region_bbox = active_region_bbox
        self.quantum_dot_polygons = quantum_dot_polygons
        self.layer = layer
        self.route_requests = []
        self.routes = []

        # 线宽定义 (nm -> um)
        self.W = {
            'tier1': 0.040,  # 40nm
            'tier2': 0.100,  # 100nm
            'tier3': 0.500   # 500nm
        }
        # 间距定义 (间距必须大于线宽)
        # 假设间距至少是线宽的1.5倍，或者用户指定的最小间距
        self.S = {
            'tier1': 0.060,  # 40nm线之间间距 > 40nm
            'tier2': 0.150,  # 100nm线之间间距 > 100nm
            'tier3': 0.750   # 500nm线之间间距 > 500nm
        }
        # 初始引出线长度和pad连接线长度
        self.initial_extension_len = 0.500  # 500nm长度的40nm宽度线
        self.pad_connection_extension_len = pad_connection_extension_len  # 500nm线连接pad引线的部分，垂直引出10um (default=5um, user specified 10um)

        # S-curve parameters (removed for now, focusing on orthogonal)
        # self.s_curve_initial_straight_len = 0.5
        # self.s_curve_45_len = 0.5
        # self.s_curve_final_straight_len = 0.5

        # 用于布线规划的参数 (可以根据需要调整)
        # Removed self.tier_radial_distances as per user request
        self.min_segment_length = 0.01 # 最小线段长度，避免生成过短的线段

    def add_route_request(self, label, start_pt, end_pt, side_info, electrode_direction=None):
        self.route_requests.append({
            'label': label,
            'start_pt': start_pt,
            'end_pt': end_pt,
            'side_info': side_info,
            'electrode_direction': electrode_direction
        })

    def plan_and_generate_routes(self):
        print("\n--- 开始全局布线规划 ---")
        all_geometries = []
        
        # Group requests by side for organized routing
        requests_by_side = {'top': [], 'bottom': [], 'left': [], 'right': []}
        for req in self.route_requests:
            requests_by_side[req['side_info'][0]].append(req)

        for side, requests in requests_by_side.items():
            if not requests:
                continue
            # Sort requests to ensure consistent routing order
            if side in ['top', 'bottom']:
                requests.sort(key=lambda r: r['start_pt'][0]) # Sort by X for top/bottom
            else:
                requests.sort(key=lambda r: r['start_pt'][1], reverse=True) # Sort by Y for left/right

            for i, req in enumerate(requests):
                label = req['label']
                start_pt = np.array(req['start_pt'])
                end_pt = np.array(req['end_pt'])
                electrode_direction = req['electrode_direction']

                # Calculate target waypoints for the current route
                target_waypoints = self._calculate_target_waypoints(start_pt, end_pt, side, i, len(requests), electrode_direction)

                # Generate path segments with varying widths
                full_path_points = []
                full_path_widths = []

                # Segment 1: Initial extension (40nm) - now orthogonal
                segment1_points = self._create_orthogonal_or_45_path(target_waypoints['p0'], target_waypoints['p0_ext'], side)
                full_path_points.extend(segment1_points)
                full_path_widths.extend([self.W['tier1']] * len(segment1_points))

                # Segment 2: 40nm to 100nm transition
                segment2_start_point = target_waypoints['p0_ext_to_p1']
                full_path_points.append(segment2_start_point)
                full_path_widths.append(self.W['tier1']) # End of 40nm segment
                full_path_points.append(segment2_start_point)
                full_path_widths.append(self.W['tier2']) # Start of 100nm segment

                segment2_points = self._create_orthogonal_or_45_path(segment2_start_point, target_waypoints['p1_to_p2'], side)
                full_path_points.extend(segment2_points[1:]) # Exclude the first point as it's already added
                full_path_widths.extend([self.W['tier2']] * (len(segment2_points) - 1))

                # Segment 3: 100nm to 500nm transition
                segment3_start_point = target_waypoints['p1_to_p2']
                full_path_points.append(segment3_start_point)
                full_path_widths.append(self.W['tier2']) # End of 100nm segment
                full_path_points.append(segment3_start_point)
                full_path_widths.append(self.W['tier3']) # Start of 500nm segment

                segment3_points = self._create_orthogonal_or_45_path(segment3_start_point, target_waypoints['p2_to_p3'], side)
                full_path_points.extend(segment3_points[1:])
                full_path_widths.extend([self.W['tier3']] * (len(segment3_points) - 1))

                # Segment 4: 500nm to pad connection
                segment4_points = self._create_orthogonal_or_45_path(target_waypoints['p2_to_p3'], target_waypoints['p3_ext'], side)
                full_path_points.extend(segment4_points[1:])
                full_path_widths.extend([self.W['tier3']] * (len(segment4_points) - 1))

                # Final connection to pad
                segment_final_points = self._create_orthogonal_or_45_path(target_waypoints['p3_ext'], target_waypoints['p_final'], side)
                full_path_points.extend(segment_final_points[1:])
                full_path_widths.extend([self.W['tier3']] * (len(segment_final_points) - 1))

                # Create FlexPath
                if len(full_path_points) > 1 and len(full_path_points) == len(full_path_widths):
                    full_route_path = gdstk.FlexPath(full_path_points, full_path_widths, joins="miter", ends="flush", layer=self.layer)
                    merged_route_geo = gdstk.boolean(full_route_path.to_polygons(), [], 'or', layer=self.layer)
                    if merged_route_geo:
                        # Store the maximum width used in this route for accurate DRC spacing checks
                        max_width_in_route = max(full_path_widths)
                        self.routes.append({'geo': merged_route_geo[0], 'label': label, 'max_width': max_width_in_route})
                        all_geometries.extend(merged_route_geo)
                else:
                    print(f"Warning: Skipping route {label} due to insufficient points or width mismatch.")

        print("全局规划和几何生成完毕。")
        return all_geometries

    def _calculate_target_waypoints(self, p0, p_final, side, index, total_on_side, electrode_direction):
        waypoints = {'p0': p0, 'p_final': p_final}
        active_min_x, active_min_y = self.active_region_bbox[0]
        active_max_x, active_max_y = self.active_region_bbox[1]

        # Calculate available space for routing on this side
        if side in ['top', 'bottom']:
            available_space = active_max_x - active_min_x
            active_start_coord = active_min_x
            p0_coord = p0[0]
        else: # 'left', 'right'
            available_space = active_max_y - active_min_y
            active_start_coord = active_min_y
            p0_coord = p0[1]

        # Calculate lateral offset based on available space and number of routes
        min_route_spacing = self.W['tier3'] + self.S['tier3'] # Use largest tier for minimum spacing
        
        # Distribute routes evenly across the available space
        if total_on_side > 1:
            total_width_needed = total_on_side * self.W['tier3'] + (total_on_side - 1) * min_route_spacing
            if total_width_needed < available_space:
                lateral_position_in_active = active_start_coord + (available_space - total_width_needed) / 2.0 + index * (self.W['tier3'] + min_route_spacing) + self.W['tier3'] / 2.0
            else:
                lateral_position_in_active = active_start_coord + (index + 0.5) * (available_space / total_on_side)
        else:
            lateral_position_in_active = active_start_coord + available_space / 2.0

        # 1. p0_initial_ext: Initial 0.5um, 40nm wide extension
        p0_initial_ext = np.array(p0)
        if electrode_direction == 'down':
            p0_initial_ext[1] -= self.initial_extension_len
        elif electrode_direction == 'up':
            p0_initial_ext[1] += self.initial_extension_len
        elif electrode_direction == 'left':
            p0_initial_ext[0] -= self.initial_extension_len
        elif electrode_direction == 'right':
            p0_initial_ext[0] += self.initial_extension_len
        else: # Default extension based on side if no specific electrode direction
            if side == 'top': p0_initial_ext[1] += self.initial_extension_len
            elif side == 'bottom': p0_initial_ext[1] -= self.initial_extension_len
            elif side == 'left': p0_initial_ext[0] -= self.initial_extension_len
            elif side == 'right': p0_initial_ext[0] += self.initial_extension_len
        
        if side in ['top', 'bottom']:
            p0_initial_ext[0] = lateral_position_in_active
        else:
            p0_initial_ext[1] = lateral_position_in_active
        waypoints['p0_ext'] = p0_initial_ext

        # Define fixed target coordinates for each tier to ensure separation and avoid overlaps
        # These values are chosen to be well within the active region and provide enough space
        # for transitions and to avoid overlaps.
        if side == 'top':
            tier1_target_y = active_min_y + (active_max_y - active_min_y) * 0.25 # 25% into active region
            tier2_target_y = active_min_y + (active_max_y - active_min_y) * 0.50 # 50% into active region
            tier3_target_y = active_min_y + (active_max_y - active_min_y) * 0.75 # 75% into active region
        elif side == 'bottom':
            tier1_target_y = active_max_y - (active_max_y - active_min_y) * 0.25
            tier2_target_y = active_max_y - (active_max_y - active_min_y) * 0.50
            tier3_target_y = active_max_y - (active_max_y - active_min_y) * 0.75
        elif side == 'left':
            tier1_target_x = active_max_x - (active_max_x - active_min_x) * 0.25
            tier2_target_x = active_max_x - (active_max_x - active_min_x) * 0.50
            tier3_target_x = active_max_x - (active_max_x - active_min_x) * 0.75
        elif side == 'right':
            tier1_target_x = active_min_x + (active_max_x - active_min_x) * 0.25
            tier2_target_x = active_min_x + (active_max_x - active_min_x) * 0.50
            tier3_target_x = active_min_x + (active_max_x - active_min_x) * 0.75

        # 2. p1_tier1: First tier waypoint (40nm)
        p1_tier1 = np.array(p0_initial_ext)
        if side in ['top', 'bottom']:
            p1_tier1[1] = tier1_target_y
            p1_tier1[0] = lateral_position_in_active
        else: # 'left', 'right'
            p1_tier1[0] = tier1_target_x
            p1_tier1[1] = lateral_position_in_active
        waypoints['p0_ext_to_p1'] = p1_tier1

        # 3. p2_tier2: Second tier waypoint (100nm)
        p2_tier2 = np.array(p1_tier1)
        if side in ['top', 'bottom']:
            p2_tier2[1] = tier2_target_y
            p2_tier2[0] = lateral_position_in_active
        else: # 'left', 'right'
            p2_tier2[0] = tier2_target_x
            p2_tier2[1] = lateral_position_in_active
        waypoints['p1_to_p2'] = p2_tier2

        # 4. p3_tier3: Third tier waypoint (500nm)
        p3_tier3 = np.array(p2_tier2)
        if side in ['top', 'bottom']:
            p3_tier3[1] = tier3_target_y
            p3_tier3[0] = lateral_position_in_active
        else: # 'left', 'right'
            p3_tier3[0] = tier3_target_x
            p3_tier3[1] = lateral_position_in_active
        waypoints['p2_to_p3'] = p3_tier3

        # 5. p_pad_ext: Orthogonal extension from p3_tier3 towards the pad
        p_pad_ext = np.array(p_final)
        if side == 'top':
            p_pad_ext[1] -= self.pad_connection_extension_len
        elif side == 'bottom':
            p_pad_ext[1] += self.pad_connection_extension_len
        elif side == 'left':
            p_pad_ext[0] += self.pad_connection_extension_len
        elif side == 'right':
            p_pad_ext[0] -= self.pad_connection_extension_len
        waypoints['p3_ext'] = p_pad_ext

        return waypoints

    def _create_orthogonal_or_45_path(self, p_start, p_end, side, is_s_curve_segment=False):
        # This function generates a path segment between p_start and p_end
        # ensuring orthogonal or 45-degree turns.
        # It returns a list of points including p_start and p_end, and any intermediate points.

        path_segment = [p_start]
        
        dx = p_end[0] - p_start[0]
        dy = p_end[1] - p_start[1]

        # If points are identical, return
        if np.isclose(dx, 0) and np.isclose(dy, 0):
            return path_segment
        
        # If it's an S-curve segment and not already orthogonal/45-degree
        if is_s_curve_segment and not (np.isclose(dx, 0) or np.isclose(dy, 0) or np.isclose(abs(dx), abs(dy))):
            # Simplified S-curve with two 45-degree turns
            s_initial_straight = self.s_curve_initial_straight_len
            s_45_len = self.s_curve_45_len
            s_final_straight = self.s_curve_final_straight_len
            
            if side == 'top':
                # Initial straight segment
                p_mid1 = np.array([p_start[0], p_start[1] + s_initial_straight])
                # First 45-degree turn
                h_dir = 1 if p_end[0] > p_mid1[0] else -1
                p_mid2 = np.array([p_mid1[0] + h_dir * s_45_len, p_mid1[1] + h_dir * s_45_len])
                # Second 45-degree turn
                p_mid3 = np.array([p_mid2[0] + h_dir * s_45_len, p_mid2[1] - h_dir * s_45_len])
                # Final straight segment
                p_mid4 = np.array([p_mid3[0], p_end[1] - s_final_straight])
                path_segment.extend([p_mid1, p_mid2, p_mid3, p_mid4, p_end])
                return path_segment
            
            elif side == 'bottom':
                # Initial straight segment
                p_mid1 = np.array([p_start[0], p_start[1] - s_initial_straight])
                # First 45-degree turn
                h_dir = 1 if p_end[0] > p_mid1[0] else -1
                p_mid2 = np.array([p_mid1[0] + h_dir * s_45_len, p_mid1[1] - h_dir * s_45_len])
                # Second 45-degree turn
                p_mid3 = np.array([p_mid2[0] + h_dir * s_45_len, p_mid2[1] + h_dir * s_45_len])
                # Final straight segment
                p_mid4 = np.array([p_mid3[0], p_end[1] + s_final_straight])
                path_segment.extend([p_mid1, p_mid2, p_mid3, p_mid4, p_end])
                return path_segment
            
            elif side == 'left':
                # Initial straight segment
                p_mid1 = np.array([p_start[0] - s_initial_straight, p_start[1]])
                # First 45-degree turn
                v_dir = 1 if p_end[1] > p_mid1[1] else -1
                p_mid2 = np.array([p_mid1[0] - v_dir * s_45_len, p_mid1[1] + v_dir * s_45_len])
                # Second 45-degree turn
                p_mid3 = np.array([p_mid2[0] + v_dir * s_45_len, p_mid2[1] + v_dir * s_45_len])
                # Final straight segment
                p_mid4 = np.array([p_end[0] + s_final_straight, p_mid3[1]])
                path_segment.extend([p_mid1, p_mid2, p_mid3, p_mid4, p_end])
                return path_segment
            
            elif side == 'right':
                # Initial straight segment
                p_mid1 = np.array([p_start[0] + s_initial_straight, p_start[1]])
                # First 45-degree turn
                v_dir = 1 if p_end[1] > p_mid1[1] else -1
                p_mid2 = np.array([p_mid1[0] + v_dir * s_45_len, p_mid1[1] + v_dir * s_45_len])
                # Second 45-degree turn
                p_mid3 = np.array([p_mid2[0] - v_dir * s_45_len, p_mid2[1] + v_dir * s_45_len])
                # Final straight segment
                p_mid4 = np.array([p_end[0] - s_final_straight, p_mid3[1]])
                path_segment.extend([p_mid1, p_mid2, p_mid3, p_mid4, p_end])
                return path_segment
            
        # Orthogonal or 45-degree direct path (L-bend or straight)
        if np.isclose(dx, 0) and np.isclose(dy, 0): # Points are identical
            return path_segment
        elif np.isclose(dx, 0) or np.isclose(dy, 0): # Orthogonal
            path_segment.append(p_end)
            return path_segment
        
        # Check for 45-degree alignment
        if np.isclose(abs(dx), abs(dy)):
            path_segment.append(p_end)
            return path_segment

        # Try to create an L-bend (orthogonal path)
        # Prioritize direction based on the side of the chip the route is on
        if side in ['top', 'bottom']: # Primarily vertical routing
            # First move horizontally, then vertically
            mid_pt = (p_end[0], p_start[1])
            path_segment.extend([mid_pt, p_end])
        else: # Primarily horizontal routing (left, right)
            # First move vertically, then horizontally
            mid_pt = (p_start[0], p_end[1])
            path_segment.extend([mid_pt, p_end])
        
        return path_segment

    def check_drc(self):
        print("\n--- 开始执行DRC检查 ---")
        violations_found = False
        num_routes = len(self.routes)

        if num_routes < 2:
            print("路由少于2条，无需检查。")
            return False

        # 1. Angle Check (Orthogonal or 45°)
        print("进行角度检查...")
        for route_info in self.routes:
            geo = route_info['geo']
            if isinstance(geo, gdstk.Polygon):
                vertices = geo.points
                for k in range(len(vertices)):
                    p1 = vertices[k-1]
                    p2 = vertices[k]
                    p3 = vertices[(k+1)%len(vertices)]
                    v1 = p2 - p1
                    v2 = p3 - p2
                    
                    # Avoid division by zero for zero-length segments
                    if np.linalg.norm(v1) < 1e-9 or np.linalg.norm(v2) < 1e-9:
                        continue

                    angle_rad = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
                    angle_deg = np.degrees(angle_rad)

                    # Check if angle is close to 0, 45, 90, 135, 180 degrees
                    is_valid_angle = False
                    for target_angle in [0, 45, 90, 135, 180]:
                        if np.isclose(angle_deg, target_angle, atol=5) or np.isclose(angle_deg, 360 - target_angle, atol=5):
                            is_valid_angle = True
                            break
                    
                    if not is_valid_angle:
                        print(f"[DRC角度违规] 路由 '{route_info['label']}' 发现一个角度为 {angle_deg:.2f}° 的非正交或45度转角。 涉及顶点: {p1}, {p2}, {p3}")
                        violations_found = True
            else:
                print(f"Warning: Route '{route_info['label']}' geometry is not a Polygon, skipping angle check.")

        # 2. Spacing Check (between different lines)
        print("进行间距检查...")
        for i in range(num_routes):
            for j in range(i + 1, num_routes):
                route_i_info = self.routes[i]
                route_j_info = self.routes[j]
                route_i_geo = route_i_info['geo']
                route_j_geo = route_j_info['geo']

                width_i = self._get_route_max_width(route_i_info)
                width_j = self._get_route_max_width(route_j_info)
                min_spacing = max(width_i, width_j) # Spacing must be greater than line width

                # Create a buffer around one route and check for intersection with the other
                # The buffer should be min_spacing / 2.0 on each side to ensure total spacing
                buffered_route_i = gdstk.offset([route_i_geo], min_spacing / 2.0, join="miter", layer=self.layer + 1)
                
                if buffered_route_i:
                    overlap_check = gdstk.boolean(buffered_route_i, [route_j_geo], 'and', layer=self.layer + 2)
                    if overlap_check:
                        violations_found = True
                        overlap_bbox = overlap_check[0].bounding_box() if overlap_check else "N/A"
                        print(f"[DRC间距违规] 路由 '{route_i_info['label']}' (宽度: {width_i:.3f}um) 和 '{route_j_info['label']}' (宽度: {width_j:.3f}um) 之间间距不足或重叠。 重叠区域: {overlap_bbox}")

        # 3. Active Region and Quantum Dot Overlap Check
        print("进行区域和量子点重叠检查...")
        active_region_poly = gdstk.rectangle(self.active_region_bbox[0], self.active_region_bbox[1], layer=self.layer + 3)
        
        for route_info in self.routes:
            route_geo = route_info['geo']
            
            # Check if route is within active region
            route_outside_active = gdstk.boolean([route_geo], [active_region_poly], 'not', layer=self.layer + 4)
            if route_outside_active:
                violations_found = True
                outside_bbox = route_info['geo'].bounding_box() if route_info['geo'] else "N/A"
                print(f"[DRC区域违规] 路由 '{route_info['label']}' 未完全在Active Region内。 外部区域: {outside_bbox}")
            
            # Check for overlap with quantum dots
            for qd_poly in self.quantum_dot_polygons:
                overlap_qd = gdstk.boolean([route_info['geo']], [qd_poly], 'and', layer=self.layer + 5)
                if overlap_qd:
                    violations_found = True
                    overlap_bbox = overlap_qd[0].bounding_box() if overlap_qd else "N/A"
                    # Attempt to find the label of the overlapping quantum dot if available
                    qd_label = "Unknown QD"
                    # This would require passing QD labels to the Router, which is not currently done.
                    # For now, we'll just report the bounding box of the QD.
                    qd_bbox = qd_poly.bounding_box()
                    print(f"[DRC重叠违规] 路由 '{route_info['label']}' 与量子点 (bbox: {qd_bbox}) 重叠。 重叠区域: {overlap_bbox}")

        if not violations_found:
            print("[DRC通过] 未发现违规。")
        return violations_found

    def _get_route_max_width(self, route_info):
        # Retrieve the stored maximum width for the route
        if 'max_width' in route_info:
            return route_info['max_width']
        
        # Fallback to bounding box heuristic if max_width is not stored (should not happen with new logic)
        if isinstance(route_info['geo'], gdstk.Polygon):
            bbox = route_info['geo'].bounding_box()
            if bbox:
                (min_x, min_y), (max_x, max_y) = bbox
                # For a polygon, a better estimate of width might be needed,
                # but for now, we rely on 'max_width' being stored.
                # If not available, this fallback is a rough estimate.
                return min(max_x - min_x, max_y - min_y)
        return 0.0 # Default to 0 if not a polygon or no bbox

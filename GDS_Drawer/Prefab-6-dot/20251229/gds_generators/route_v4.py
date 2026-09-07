import gdstk
import numpy as np
import math
import re

class Route_v4:
    """
    一个代表单条完整扇出引线的类，它接收预先计算好的、DRC-clean的路径点，
    并严格按照手绘草图的样式构建独立的几何分段。
    """
    def __init__(self, waypoints, layer):
        self.waypoints = waypoints
        self.layer = layer
        # Updated wire widths: 40nm (0.04), 100nm (0.1), 1um (1.0)
        self.W = {1: 0.04, 2: 0.1, 3: 1.0} 
        self.geometries = []

    def _build_geometry(self):
        """根据预计算的waypoint，分段构建独立的几何对象。"""
        # --- 段 1 (40nm): 正交引出 ---
        p0, p1 = self.waypoints['p0'], self.waypoints['p1']
        self.geometries.append(gdstk.FlexPath([p0, p1], self.W[1], joins="miter", ends="flush"))
        
        # --- 段 2 (100nm): 带正交桩的45度S型走线 ---
        p2, p1_stub_end, m1, m2, p2_stub_start = self.waypoints['p2'], self.waypoints['p1_stub_end'], self.waypoints['m1'], self.waypoints['m2'], self.waypoints['p2_stub_start']
        self.geometries.append(gdstk.FlexPath([p1, p1_stub_end], self.W[2], joins="miter", ends="flush"))
        self.geometries.append(gdstk.FlexPath([p1_stub_end, m1, m2, p2_stub_start], self.W[2], joins="miter", ends="flush"))
        self.geometries.append(gdstk.FlexPath([p2_stub_start, p2], self.W[2], joins="miter", ends="flush"))

        # --- 段 3 (1um): 带正交桩的45度S型走线 ---
        p3, p2_stub_end, m3, m4, p3_stub_start = self.waypoints['p3'], self.waypoints['p2_stub_end'], self.waypoints['m3'], self.waypoints['m4'], self.waypoints['p3_stub_start']
        self.geometries.append(gdstk.FlexPath([p2, p2_stub_end], self.W[3], joins="miter", ends="flush"))
        self.geometries.append(gdstk.FlexPath([p2_stub_end, m3, m4, p3_stub_start], self.W[3], joins="miter", ends="flush"))
        self.geometries.append(gdstk.FlexPath([p3_stub_start, p3], self.W[3], joins="miter", ends="flush"))
        
        # --- 段 4 (1um): 最终垂直对接 ---
        p_final = self.waypoints['p_final']
        if abs(p3[0] - p_final[0]) < abs(p3[1] - p_final[1]): # Vertical final segment
            self.geometries.append(gdstk.FlexPath([p3, (p3[0], p_final[1]), p_final], self.W[3], joins="miter", ends="flush"))
        else: # Horizontal final segment
            self.geometries.append(gdstk.FlexPath([p3, (p_final[0], p3[1]), p_final], self.W[3], joins="miter", ends="flush"))

    def get_geometry(self):
        """计算并返回最终合并后的几何形状。"""
        self._build_geometry()
        all_polys = [poly for path in self.geometries for poly in path.to_polygons()]
        return gdstk.boolean(all_polys, [], 'or', layer=self.layer)

class RouteManager_v4:
    """
    一个高级的布线管理器，负责全局规划和DRC检查。
    """
    def __init__(self, chip_center, layer):
        self.chip_center = chip_center
        self.layer = layer
        self.route_requests = []
        self.routes = []

        # Updated wire widths: 40nm (0.04), 100nm (0.1), 1um (1.0)
        self.W = {1: 0.04, 2: 0.1, 3: 1.0}
        self.S = {1: 0.08, 2: 0.2, 3: 2.0} # Increased spacing to mitigate overlap issues
        self.P = {1: self.W[1]+self.S[1], 2: self.W[2]+self.S[2], 3: self.W[3]+self.S[3]}
        self.D = {1: 15, 2: 40, 3: 70}
        self.stub_len = {1: 5.0, 2: 10.0} # Increased stub length for smoother S-bends

    def add_route_request(self, label, start_pt, end_pt, side_info):
        self.route_requests.append({'label': label, 'start_pt': start_pt, 'end_pt': end_pt, 'side_info': side_info})

    def plan_and_generate_routes(self):
        print("\n--- 开始全局布线规划 ---")
        requests_by_side = {'top': [], 'bottom': [], 'left': [], 'right': []}
        for req in self.route_requests: requests_by_side[req['side_info'][0]].append(req)

        all_geometries = []
        for side, requests in requests_by_side.items():
            if not requests: continue
            requests.sort(key=lambda r: r['side_info'][1])
            total = len(requests)
            for i, req in enumerate(requests):
                req['side_info'] = (side, i, total)
                waypoints = self._calculate_waypoints_for_route(req)
                router = Route_v4(waypoints=waypoints, layer=self.layer)
                geometries = router.get_geometry()
                if geometries:
                    self.routes.append({'geo': geometries[0], 'label': req.get('label', f'route_{len(self.routes)}')})
                all_geometries.extend(geometries)
        
        print("全局规划和几何生成完毕。")
        return all_geometries

    def _calculate_waypoints_for_route(self, req):
        waypoints = {'p0': req['start_pt'], 'p_final': req['end_pt']}
        side, index, total = req['side_info']
        center = self.chip_center

        for tier in [1, 2, 3]:
            pitch = self.P[tier]
            offset = (index - (total - 1) / 2.0) * pitch
            if side == 'bottom':   waypoints[f'p{tier}'] = np.array((center[0] + offset, center[1] - self.D[tier]))
            elif side == 'top':    waypoints[f'p{tier}'] = np.array((center[0] + offset, center[1] + self.D[tier]))
            elif side == 'left':   waypoints[f'p{tier}'] = np.array((center[0] - self.D[tier], center[1] - offset))
            else: # right
                waypoints[f'p{tier}'] = np.array((center[0] + self.D[tier], center[1] - offset))

        p1, p2, p3 = waypoints['p1'], waypoints['p2'], waypoints['p3']
        
        waypoints['p1_stub_end'], waypoints['m1'], waypoints['m2'], waypoints['p2_stub_start'] = self._create_s_bend_points(p1, p2, side, self.stub_len[1])
        waypoints['p2_stub_end'], waypoints['m3'], waypoints['m4'], waypoints['p3_stub_start'] = self._create_s_bend_points(p2, p3, side, self.stub_len[2])
        return waypoints
    
    def _create_s_bend_points(self, p_start, p_end, side, stub_len):
        p_start = np.array(p_start)
        p_end = np.array(p_end)

        dx_total = p_end[0] - p_start[0]
        dy_total = p_end[1] - p_start[1]

        # Determine the initial and final orthogonal segments
        if side in ['top', 'bottom']: # Main routing direction is horizontal, S-bend is vertical
            p_start_ortho = p_start + np.array([0, np.sign(dy_total) * stub_len])
            p_end_ortho = p_end - np.array([0, np.sign(dy_total) * stub_len])
        else: # Main routing direction is vertical, S-bend is horizontal
            p_start_ortho = p_start + np.array([np.sign(dx_total) * stub_len, 0])
            p_end_ortho = p_end - np.array([np.sign(dx_total) * stub_len, 0])

        # Calculate intermediate points for the 45-degree S-bend
        # The S-bend will consist of two 45-degree segments and one orthogonal segment in between.
        # The length of the diagonal part will be `stub_len`

        # Calculate intermediate points for the 45-degree S-bend
        # The S-bend will consist of two 45-degree segments and one orthogonal segment in between.
        # The `stub_len` defines the length of the orthogonal projection for the 45-degree turns.

        s_dx = np.sign(dx_total)
        s_dy = np.sign(dy_total)

        if side in ['top', 'bottom']: # Main routing direction is horizontal, S-bend is vertical
            # The middle segment will be horizontal.
            # Calculate the y-coordinate for the horizontal middle segment.
            # This y-coordinate should be between p_start_ortho[1] and p_end_ortho[1].
            # The 45-degree turns will extend 'stub_len' in x and y.
            
            # Determine the y-coordinate of the first 45-degree turn's end point (m1_y)
            # and the second 45-degree turn's start point (m2_y).
            # These should be consistent to form a horizontal middle segment.
            
            # The available vertical distance for the 45-degree turns.
            available_dy = abs(p_end_ortho[1] - p_start_ortho[1])
            
            # The actual length of the diagonal projection, limited by available space.
            # It should be at most half of the available vertical distance.
            diagonal_len = min(stub_len, available_dy / 2.0)
            
            m1_y = p_start_ortho[1] + s_dy * diagonal_len
            m2_y = p_end_ortho[1] - s_dy * diagonal_len
            
            # If the diagonal_len is too small, the middle segment might not be long enough.
            # Ensure m1_y and m2_y are the same for a horizontal middle segment.
            # This means the actual diagonal_len used for the y-component is (p_end_ortho[1] - p_start_ortho[1]) / (2 * s_dy)
            # Let's use a common y for the horizontal segment.
            common_y = (p_start_ortho[1] + p_end_ortho[1]) / 2 # Midpoint of the vertical span
            
            m1 = np.array([p_start_ortho[0] + s_dx * diagonal_len, common_y])
            m2 = np.array([p_end_ortho[0] - s_dx * diagonal_len, common_y])

            # Ensure m1_x and m2_x don't cross. If they do, adjust the x-coordinates to meet at the midpoint.
            if s_dx * (m2[0] - m1[0]) < 0:
                mid_x_segment = (m1[0] + m2[0]) / 2
                m1[0] = mid_x_segment - s_dx * diagonal_len
                m2[0] = mid_x_segment + s_dx * diagonal_len

        else: # Main routing direction is vertical, S-bend is horizontal
            # The middle segment will be vertical.
            available_dx = abs(p_end_ortho[0] - p_start_ortho[0])
            diagonal_len = min(stub_len, available_dx / 2.0)

            m1_x = p_start_ortho[0] + s_dx * diagonal_len
            m2_x = p_end_ortho[0] - s_dx * diagonal_len

            common_x = (p_start_ortho[0] + p_end_ortho[0]) / 2 # Midpoint of the horizontal span

            m1 = np.array([common_x, p_start_ortho[1] + s_dy * diagonal_len])
            m2 = np.array([common_x, p_end_ortho[1] - s_dy * diagonal_len])

            # Ensure m1_y and m2_y don't cross. If they do, adjust the y-coordinates to meet at the midpoint.
            if s_dy * (m2[1] - m1[1]) < 0:
                mid_y_segment = (m1[1] + m2[1]) / 2
                m1[1] = mid_y_segment - s_dy * diagonal_len
                m2[1] = mid_y_segment + s_dy * diagonal_len
        
        return p_start_ortho, m1, m2, p_end_ortho

    def check_drc(self):
        print("\n--- 开始执行DRC检查 ---"); violations_found = False; num_routes = len(self.routes)
        if num_routes < 2: print("路由少于2条，无需检查。"); return False
        
        # Angle Check
        all_polys = [r['geo'] for r in self.routes]
        merged = gdstk.boolean(all_polys, [], 'or')
        for poly in merged:
            vertices = poly.points
            for k in range(len(vertices)):
                p1 = vertices[k-1]; p2 = vertices[k]; p3 = vertices[(k+1)%len(vertices)]
                v1 = p2 - p1; v2 = p3 - p2
                if np.linalg.norm(v1) * np.linalg.norm(v2) > 1e-12:
                    angle = np.degrees(np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))
                    # Check for angles that are not multiples of 45 degrees (with a tolerance)
                    if not np.isclose(angle % 45, 0, atol=1) and not np.isclose(angle % 45, 45, atol=1):
                        print(f"[DRC角度违规] 发现一个角度为 {angle:.2f}° 的非45度倍数转角。")
                        violations_found = True
        
        # Overlap Check
        for i in range(num_routes):
            for j in range(i + 1, num_routes):
                overlap = gdstk.boolean(self.routes[i]['geo'], self.routes[j]['geo'], 'and', layer=self.layer + 10)
                if overlap: 
                    violations_found = True
                    for poly in overlap:
                        (min_x, min_y), (max_x, max_y) = poly.bounding_box()
                        print(f"[DRC重叠违规] 在路由 {i} ({self.routes[i]['label']}) 和 {j} ({self.routes[j]['label']}) 之间发现重叠! 坐标: ({(min_x+max_x)/2:.2f}, {(min_y+max_y)/2:.2f})")

        if not violations_found: print("[DRC通过] 未发现违规。")
        return violations_found

import gdsfactory as gf
from gdsfactory.component import Component
from gdsfactory.port import Port
from typing import List, Tuple, Optional
from gdsfactory.routing import route_quad
from gdsfactory.routing.route_single import route_single

# Removed @cell decorator to resolve 'TypeError: unhashable type: 'DPort''
def route_v3(
    quantum_dot_ports: List[Port],
    pad_ports: List[Port],
    layer: int = 1,
    width_40nm: float = 0.04,  # 40nm in um
    width_100nm: float = 0.1,  # 100nm in um
    width_500nm: float = 0.5,  # 500nm in um
    quantum_dot_extension: float = 1.0, # Distance to extend from QD before S-bend
    s_bend_radius: float = 5.0, # Radius for S-bends
    min_100nm_spacing: float = 0.1, # Minimum spacing for 100nm lines
    min_500nm_spacing: float = 0.5, # Minimum spacing for 500nm lines
) -> Component:
    """
    Routes from quantum dot electrodes to pad leads with varying widths.

    Args:
        quantum_dot_ports: List of ports from quantum dot electrodes.
        pad_ports: List of ports from pad leads.
        layer: GDS layer for routing.
        width_40nm: Width of the quantum dot electrode connection (in um).
        width_100nm: Width of the intermediate connection (in um).
        width_500nm: Width of the pad lead connection (in um).
        quantum_dot_extension: Distance to extend horizontally/vertically from QD before S-bend.
        s_bend_radius: Radius for S-bends.
        min_100nm_spacing: Minimum spacing requirement for 100nm lines.
        min_500nm_spacing: Minimum spacing requirement for 500nm lines.

    Returns:
        A gdsfactory Component with the routed connections.
    """
    c = Component()

    # Ensure the number of quantum dot ports matches pad ports for 1:1 routing
    if len(quantum_dot_ports) != len(pad_ports):
        raise ValueError("Number of quantum dot ports must match number of pad ports.")

    # Sort ports to ensure consistent routing order (e.g., by y-coordinate)
    quantum_dot_ports.sort(key=lambda p: p.center[1])
    pad_ports.sort(key=lambda p: p.center[1])

    for i, (qd_port, pad_port) in enumerate(zip(quantum_dot_ports, pad_ports)):
        # Current port for routing
        current_port = qd_port

        # 1. 40nm section: Straight extension
        # Calculate the end point of the initial straight extension
        if current_port.orientation == 0: # East
            end_40nm_straight = (current_port.center[0] + quantum_dot_extension, current_port.center[1])
        elif current_port.orientation == 180: # West
            end_40nm_straight = (current_port.center[0] - quantum_dot_extension, current_port.center[1])
        elif current_port.orientation == 90: # North
            end_40nm_straight = (current_port.center[0], current_port.center[1] + quantum_dot_extension)
        elif current_port.orientation == 270: # South
            end_40nm_straight = (current_port.center[0], current_port.center[1] - quantum_dot_extension)
        else:
            raise ValueError(f"Unsupported quantum dot port orientation: {current_port.orientation}")

        # Add the initial straight segment
        straight_segment_comp = gf.path.Path([current_port.center, end_40nm_straight]).extrude(width=width_40nm, layer=layer)
        c.add_ref(straight_segment_comp)

        # Update current_port
        current_port = gf.Port(
            name=f"qd_ext_straight_end_{i}",
            center=end_40nm_straight,
            width=width_40nm,
            orientation=current_port.orientation,
            layer=layer
        )

        # 1.5. 40nm S-bend
        # Use gf.components.bend_s
        s_bend_dx = s_bend_radius * 4 # Length of the S-bend in the direction of travel
        s_bend_dy = s_bend_radius * 2 * (1 if i % 2 == 0 else -1) # Lateral displacement, alternating

        s_bend_40nm_comp = gf.components.bend_s(
            size=(s_bend_dx, s_bend_dy),
            width=width_40nm
        )
        s_bend_40nm_ref = c.add_ref(s_bend_40nm_comp)
        s_bend_40nm_ref.connect("o1", current_port)
        current_port = gf.Port(
            name=s_bend_40nm_ref.ports["o2"].name,
            center=s_bend_40nm_ref.ports["o2"].center,
            width=width_40nm,
            orientation=s_bend_40nm_ref.ports["o2"].orientation,
            layer=s_bend_40nm_ref.ports["o2"].layer
        )

        # 2. Taper 40nm to 100nm
        taper_40_100_comp = gf.components.taper(
            length=s_bend_radius, width1=width_40nm, width2=width_100nm, layer=layer
        )
        taper_40_100_ref = c.add_ref(taper_40_100_comp)
        taper_40_100_ref.connect("o1", current_port)
        current_port = gf.Port(
            name=taper_40_100_ref.ports["o2"].name,
            center=taper_40_100_ref.ports["o2"].center,
            width=width_100nm,
            orientation=taper_40_100_ref.ports["o2"].orientation,
            layer=taper_40_100_ref.ports["o2"].layer
        )

        # 3. 100nm 45-degree S-bend
        # Create a path for the 100nm 45-degree S-bend
        # The S-bend needs to provide enough lateral displacement to meet min_100nm_spacing
        # Let's assume a fixed S-bend length for now, and ensure it's wide enough.
        # The lateral displacement should be at least min_100nm_spacing + width_100nm
        # For simplicity, let's use a fixed S-bend size for now.
        s_bend_100nm_dx = s_bend_radius * 4 # Length of the S-bend
        s_bend_100nm_dy = s_bend_radius * 2 * (1 if i % 2 == 0 else -1) # Lateral displacement

        s_bend_100nm_comp = gf.components.bend_s(
            size=(s_bend_100nm_dx, s_bend_100nm_dy),
            width=width_100nm
        )
        s_bend_100nm_ref = c.add_ref(s_bend_100nm_comp)
        s_bend_100nm_ref.connect("o1", current_port)
        current_port = gf.Port(
            name=s_bend_100nm_ref.ports["o2"].name,
            center=s_bend_100nm_ref.ports["o2"].center,
            width=width_100nm,
            orientation=s_bend_100nm_ref.ports["o2"].orientation,
            layer=s_bend_100nm_ref.ports["o2"].layer
        )

        # 4. Taper 100nm to 500nm
        taper_100_500_comp = gf.components.taper(
            length=s_bend_radius, width1=width_100nm, width2=width_500nm, layer=layer
        )
        taper_100_500_ref = c.add_ref(taper_100_500_comp)
        taper_100_500_ref.connect("o1", current_port)
        current_port = gf.Port(
            name=taper_100_500_ref.ports["o2"].name,
            center=taper_100_500_ref.ports["o2"].center,
            width=width_500nm,
            orientation=taper_100_500_ref.ports["o2"].orientation,
            layer=taper_100_500_ref.ports["o2"].layer
        )

        # 5. 500nm section: Orthogonal routing to pad
        # Use gf.routing.get_route for this final segment to ensure perpendicular connection.
        xs_500nm = gf.cross_section.cross_section(width=width_500nm, layer=layer, radius=s_bend_radius)

        # Define waypoints to ensure perpendicular connection.
        # The routing should end perpendicularly to the pad port.
        # This means the last segment of the route should be aligned with the pad port's orientation.
        waypoints_500nm = []
        if pad_port.orientation == 0 or pad_port.orientation == 180: # Pad is horizontal (East/West), connect vertically
            # Route to an x-coordinate that aligns with the pad, then turn to pad's y
            waypoints_500nm.append((pad_port.center[0], current_port.center[1]))
        elif pad_port.orientation == 90 or pad_port.orientation == 270: # Pad is vertical (North/South), connect horizontally
            # Route to a y-coordinate that aligns with the pad, then turn to pad's x
            waypoints_500nm.append((current_port.center[0], pad_port.center[1]))

        route_500nm_to_pad = route_single(
            port1=current_port,
            port2=pad_port,
            cross_section=xs_500nm,
            bend=gf.components.bend_circular(radius=s_bend_radius), # Use 90-degree bends
            straight=gf.components.straight(layer=layer),
        )
        # Add references and labels individually
        for ref in route_500nm_to_pad.references:
            c.add_ref(ref)
        for label in route_500nm_to_pad.labels:
            c.add(label)
    return c

if __name__ == "__main__":
    # Example usage:
    # Define a layer for the example
    example_layer = 1

    # Create dummy quantum dot ports and pad ports for testing
    qd_ports = [
        gf.Port(name="qd_e1", center=(0, 0), width=0.04, orientation=0, layer=example_layer),
        gf.Port(name="qd_e2", center=(0, 10), width=0.04, orientation=0, layer=example_layer),
    ]

    pad_ports = [
        gf.Port(name="pad_p1", center=(100, 0), width=0.5, orientation=180, layer=example_layer),
        gf.Port(name="pad_p2", center=(100, 10), width=0.5, orientation=180, layer=example_layer),
    ]

    # Generate the routing component
    routed_component = route_v3(qd_ports, pad_ports, layer=example_layer)

    # Save the component to a GDS file
    routed_component.write_gds("route_v3_project/routed_example.gds")
    print("Generated route_v3_project/routed_example.gds")

    # You can also view it in KLayout if you have it installed
    # routed_component.show()

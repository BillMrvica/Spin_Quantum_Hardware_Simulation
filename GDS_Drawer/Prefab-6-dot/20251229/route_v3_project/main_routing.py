import gdsfactory as gf
from gdsfactory.component import Component
from quantum_dot_generator import QuantumDotGenerator
from pad_lead_generator import PadLeadGenerator
from route_v3 import route_v3
from typing import Dict, List, Tuple

def main_routing_flow() -> Component:
    c = gf.Component("Full_Routed_Device")

    # 1. Generate Quantum Dot Layout
    qd_generator = QuantumDotGenerator()
    qd_component, qd_ports_dict = qd_generator.create_6qd_layout_with_labels()
    qd_ref = c.add_ref(qd_component)

    # Adjust QD ports to be relative to the main component
    qd_ports = []
    for name, port in qd_ports_dict.items():
        # Assuming qd_ref is at (0,0) for now, if placed differently, adjust port.center
        qd_ports.append(gf.Port(name=name, center=port.center, width=port.width, orientation=port.orientation, layer=port.layer))

    # 2. Generate Pad Lead Layout
    pad_generator = PadLeadGenerator()
    pad_component, pad_ports_dict, acenter, all_pads_info = pad_generator.create_rect_wire_layout()
    pad_ref = c.add_ref(pad_component)

    # Adjust pad_ref position to be far enough from QD_ref
    # For simplicity, let's place it to the right of the QD component
    qd_max_x = qd_ref.xmax
    pad_min_x = pad_ref.xmin
    x_offset = qd_max_x - pad_min_x + 200 # 200um buffer
    pad_ref.movex(x_offset)

    # Adjust pad ports to be relative to the main component after moving pad_ref
    pad_ports = []
    for name, port in pad_ports_dict.items():
        original_center = port.center
        new_center = (original_center[0] + x_offset, original_center[1])
        pad_ports.append(gf.Port(name=name, center=new_center, width=port.width, orientation=port.orientation, layer=port.layer))

    # 3. Match Quantum Dot Ports to Pad Ports
    # This is a simplified matching. In a real scenario, you'd have a more sophisticated matching logic.
    # For now, let's assume a direct mapping based on sorted y-coordinates or specific naming conventions.
    # Let's collect the relevant QD ports and Pad ports for routing.
    
    # Example: Match QD_PG1 to a pad, QD_B1 to another, etc.
    # This requires careful consideration of the naming conventions and desired connections.
    # For demonstration, let's try to match based on sorted lists.

    # Filter and sort QD ports that are likely to be routed
    qd_ports_to_route = [p for p in qd_ports if p.name.startswith('QD_') or p.name.startswith('SET')]
    qd_ports_to_route.sort(key=lambda p: p.center[1]) # Sort by y-coordinate

    # Filter and sort Pad ports that are likely to be connected
    pad_ports_to_route = [p for p in pad_ports if p.name in qd_ports_dict.keys()] # Match by name
    pad_ports_to_route.sort(key=lambda p: p.center[1]) # Sort by y-coordinate

    # Ensure equal number of ports for 1:1 routing
    if len(qd_ports_to_route) != len(pad_ports_to_route):
        print(f"Warning: Mismatch in number of QD ports ({len(qd_ports_to_route)}) and Pad ports ({len(pad_ports_to_route)}) for routing.")
        # For now, let's take the minimum number to avoid errors
        min_len = min(len(qd_ports_to_route), len(pad_ports_to_route))
        qd_ports_to_route = qd_ports_to_route[:min_len]
        pad_ports_to_route = pad_ports_to_route[:min_len]

    # 4. Generate Routing using route_v3
    if qd_ports_to_route and pad_ports_to_route:
        routed_connections = route_v3(
            quantum_dot_ports=qd_ports_to_route,
            pad_ports=pad_ports_to_route,
            layer=1, # Use layer 1 for routing
            width_40nm=0.04,
            width_100nm=0.1,
            width_500nm=0.5,
            quantum_dot_extension=5.0, # Increased extension for better visibility
            s_bend_radius=10.0, # Increased radius
        )
        c.add_ref(routed_connections)

    return c

if __name__ == "__main__":
    full_device = main_routing_flow()
    full_device.write_gds("route_v3_project/full_routed_device.gds")
    print("Generated route_v3_project/full_routed_device.gds")
    # full_device.show() # Uncomment to view in KLayout

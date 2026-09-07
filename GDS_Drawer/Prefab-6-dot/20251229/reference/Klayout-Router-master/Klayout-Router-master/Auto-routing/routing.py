import argparse
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# Ensure we can import local modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

try:
    from Krouter import Krouter
except ImportError:
    print("Error: Could not import Krouter. Make sure you are running this script from the Auto-routing directory or that Krouter.py is in the python path.")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Run Auto-Klayout routing.")
    
    # File and Cell settings
    parser.add_argument("--file_path", type=str, required=True, help="Path to the GDS file")
    parser.add_argument("--cell_name", type=str, default="Top", help="Name of the cell to route")
    
    # Layer settings
    parser.add_argument("--bbox_layer", type=str, default="7/0", help="Layer for total area bounding box")
    parser.add_argument("--pin_a_layer", type=str, default="102/1", help="Layer for Pin A")
    parser.add_argument("--pin_b_layer", type=str, default="111/1", help="Layer for Pin B")
    parser.add_argument("--output_layer", type=str, default="121/0", help="Layer to output routed paths")
    
    # Map settings
    parser.add_argument("--resolution", type=int, default=5, help="Map resolution")
    
    # Routing Configuration
    parser.add_argument("--obs_safe_distance", type=int, default=16)
    parser.add_argument("--obs_hardness", type=int, default=20)
    parser.add_argument("--obs_damping_step", type=int, default=4)
    
    parser.add_argument("--pin_a_safe_dist", type=int, default=20, dest="Pin_A_safe_distance")
    parser.add_argument("--pin_b_safe_dist", type=int, default=200, dest="Pin_B_safe_distance")
    parser.add_argument("--pin_hardness", type=int, default=20)
    parser.add_argument("--pin_a_damping", type=int, default=4, dest="Pin_A_damping_step")
    parser.add_argument("--pin_b_damping", type=int, default=1, dest="Pin_B_damping_step")
    
    parser.add_argument("--width", type=int, default=4, dest="routing_path_width", help="Routing path width")
    parser.add_argument("--path_safe_dist", type=int, default=20, dest="path_safe_distance")
    parser.add_argument("--path_hardness", type=int, default=10)
    parser.add_argument("--path_damping", type=int, default=5, dest="path_damping_step")
    parser.add_argument("--path_density_hardness", type=int, default=1)
    
    parser.add_argument("--round", type=int, default=1, help="Number of routing rounds")
    
    # Animation settings
    parser.add_argument("--animate", action="store_true", help="Enable animation generation")
    parser.add_argument("--save_anim", type=str, help="Path to save animation (e.g. routing.mp4 or routing.gif)")
    parser.add_argument("--fps", type=int, default=20, help="Frames per second for animation")

    args = parser.parse_args()

    # Obstruction rules (hardcoded default for now as per notebook)
    obs_rule1 = {"layers": ["1/0", "2/0"], "bbx": True, "pad": 6}
    
    print(f"Initializing Krouter with file: {args.file_path}")
    
    Krt = Krouter(
        file_path=args.file_path,
        cell_name=args.cell_name,
        total_area_bbox_layer=args.bbox_layer,
        Pin_A_layer=args.pin_a_layer,
        Pin_B_layer=args.pin_b_layer,
        map_resolution=args.resolution,
        obs_rules=[obs_rule1,],
    )

    print("Matching pins...")
    Krt.PinMatcher.hungarian_match()
    Krt.PinMatcher.sort_matched_pins()
    
    print("Starting routing...")
    
    # Prepare arguments for self_adaptive_route
    route_kwargs = {
        "obs_safe_distance": args.obs_safe_distance,
        "obs_hardness": args.obs_hardness,
        "obs_damping_step": args.obs_damping_step,
        "Pin_A_safe_distance": args.Pin_A_safe_distance,
        "Pin_B_safe_distance": args.Pin_B_safe_distance,
        "pin_hardness": args.pin_hardness,
        "Pin_A_damping_step": args.Pin_A_damping_step,
        "Pin_B_damping_step": args.Pin_B_damping_step,
        "routing_path_width": args.routing_path_width,
        "path_safe_distance": args.path_safe_distance,
        "path_hardness": args.path_hardness,
        "path_damping_step": args.path_damping_step,
        "path_density_hardness": args.path_density_hardness,
        "round": args.round
    }
    
    frames = Krt.self_adaptive_route(**route_kwargs)
    
    if args.animate and args.save_anim:
        print(f"Generating animation to {args.save_anim}...")
        plt.ioff() # Turn off interactive mode
        fig, ax = plt.subplots(dpi=300)
        
        # Setup plot
        im = ax.imshow(Krt.Map.map, origin='lower', cmap='YlGnBu_r', vmin=-5, vmax=30)
        fig.colorbar(im, ax=ax)
        
        def update(frame):
            im.set_data(frame)
            return [im]

        anim = FuncAnimation(fig, func=update, frames=frames, interval=200, blit=False, repeat=False)
        
        writer = "pillow"
        if args.save_anim.endswith('.mp4'):
            writer = "ffmpeg"
        
        try:
            anim.save(args.save_anim, writer=writer, fps=args.fps)
            print(f"Animation saved to {args.save_anim}")
        except Exception as e:
            print(f"Error saving animation: {e}")
            print("Make sure ffmpeg is installed if saving as mp4.")

    print(f"Inserting paths to layer {args.output_layer}...")
    Krt.insert_paths_to_layout(routing_path_layer=args.output_layer)
    
    print("Saving layout...")
    Krt.save_layout()
    print("Done.")

if __name__ == "__main__":
    main()

import schemdraw
import schemdraw.elements as elm

# Use the most basic method to create and draw the figure
d = schemdraw.Drawing(fontsize=12)

# --- Draw the main circuit component by component ---

# Input Line
d.add(elm.Line().label('(In)', loc='left'))

# First Shunt Capacitor
d.push()
d.add(elm.Capacitor().down().label('$C_1$ (fixed)'))
d.add(elm.Ground())
d.pop()

# --- First SQUID (L2) ---
# Drawn manually as two parallel JJs
start_node_L2 = d.here
d.push() # Path 1
d.add(elm.Line().up(0.25))
d.add(elm.Ic().right().label('JJ')) # Using Ic as a robust alternative for Josephson
d.add(elm.Line().down(0.25))
end_node_L2 = d.here
d.pop() # Path 2
d.add(elm.Line().at(start_node_L2).down(0.25))
d.add(elm.Ic().right().label('JJ'))
d.add(elm.Line().up(0.25).to(end_node_L2))
d.add(elm.Label().at(start_node_L2, dy=0.7).label('$L_2$ (tunable)')) # Label for the SQUID

# Second Shunt Capacitor
d.push()
d.add(elm.Capacitor().down().label('$C_3$ (fixed)'))
d.add(elm.Ground())
d.pop()

# --- Second SQUID (L4) ---
# Drawn manually again
start_node_L4 = d.here
d.push() # Path 1
d.add(elm.Line().up(0.25))
d.add(elm.Ic().right().label('JJ'))
d.add(elm.Line().down(0.25))
end_node_L4 = d.here
d.pop() # Path 2
d.add(elm.Line().at(start_node_L4).down(0.25))
d.add(elm.Ic().right().label('JJ'))
d.add(elm.Line().up(0.25).to(end_node_L4))
d.add(elm.Label().at(start_node_L4, dy=0.7).label('$L_4$ (tunable)'))

# Final Shunt Capacitor
d.push()
d.add(elm.Capacitor().down().label('$C_5$ (fixed)'))
d.add(elm.Ground())
d.pop()

# Output Line
d.add(elm.Line().right().label('(Out)', loc='right'))

# Render the final drawing
d.draw()
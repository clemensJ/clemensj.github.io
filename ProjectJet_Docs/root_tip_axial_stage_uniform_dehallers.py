#Axial Compressor Math Model with Uniform deHaller's Ratio
#This version ensures hub and tip deHaller's ratios equal the midline deHaller's ratio
import math

#INPUT PARAMETERS
RPM = 100000

s1_r = 37 # calculated meanline (mm)
inlet_area = 0.006  # m^2

drm = inlet_area*1000000/(4*math.pi*s1_r) # calculate offset from meanline to outer and inner annulus

inner_rad = s1_r-drm
outer_rad = s1_r+drm

alpha1 = 10 # IGV angle (degrees)

s1_desiredPR = 1.6

s1dh = 37000 # J/kg - specific enthalpy rise
mass_flow = 0.8 # kg/s
s1_dehallers = 0.75  # target deHaller's ratio (this will be achieved at midline, hub, and tip)

# calculated parameters
a1 = alpha1*math.pi/180
cz = mass_flow*(1/1.225)/(inlet_area) # calculate axial flow velocity through engine
Wr1 = mass_flow*s1dh # work of 1st stage (thermodynamic)

# Tangential velocities
u1 = (s1_r/1000)*(2*math.pi/60)*RPM # meanline tangential velocity
u1_hub = (inner_rad/1000)*(2*math.pi/60)*RPM # hub tangential velocity
u1_tip = (outer_rad/1000)*(2*math.pi/60)*RPM # tip tangential velocity

s1_dh_stageloading = (s1dh)/(u1*u1) # thermodynamic stage loading

# ===== STEP 1: Calculate MIDLINE angles and deHaller's ratio =====
# This will be our reference deHaller's ratio that hub and tip must match

# Calculate midline b2 from work equation
# (work rate) = mass flow * U * (C02-C01)
i = ((Wr1/(mass_flow*u1))-u1)/cz + math.tan(a1)
b2 = math.atan(i)

# Calculate midline c01 and b1
c01 = cz*math.tan(a1)
i = (u1-c01)/cz
b1 = -math.atan(i)

# Calculate midline a2
i = u1/cz + math.tan(b2)
a2 = math.atan(i)

# Calculate midline a3 (stator exit angle)
try:
    a3 = math.acos(math.cos(a2)/s1_dehallers)
except:
    a3 = 0 # purely axial
    new_dehallers = math.cos(a3)/math.cos(a2)
    print("\n*****s1 dehallers must be reduced assuming counter-swirl is undesirable: " + str(new_dehallers)+"\n\n")

c03 = cz*math.tan(a3)

# Calculate midline tangential velocity imparted
c02 = cz*math.tan(a2)

# Calculate midline enthalpy change: Δh = U * (c02 - c01)
s1dh_midline = u1 * (c02 - c01)

# Calculate midline relative velocities
w1_midline = math.sqrt(cz*cz + (u1-c01)*(u1-c01))  # W1 = rotor inlet relative speed
w2_midline = math.sqrt(cz*cz + (u1-c02)*(u1-c02))  # W2 = rotor exit relative speed

# Calculate midline deHaller's ratio (this is our target for hub and tip)
r1_dehallers_midline = w2_midline / w1_midline
print(f"Midline deHaller's ratio: {r1_dehallers_midline:.6f}")

# ===== STEP 2: Calculate HUB angles to match midline deHaller's ratio =====
# We need to solve for angles that give w2_hub/w1_hub = r1_dehallers_midline

# Hub inlet relative velocity (known from geometry)
w1_hub = math.sqrt(cz*cz + (u1_hub-c01)*(u1_hub-c01))

# For hub, we need w2_hub = r1_dehallers_midline * w1_hub
w2_hub_target = r1_dehallers_midline * w1_hub

# Now solve for c02_hub from: w2_hub² = cz² + (u1_hub - c02_hub)²
# w2_hub² - cz² = (u1_hub - c02_hub)²
# sqrt(w2_hub² - cz²) = |u1_hub - c02_hub|
# Choose the positive root (typical compressor design: c02 < u1, so u1 - c02 > 0)
w2_hub_sq = w2_hub_target * w2_hub_target
if w2_hub_sq < cz*cz:
    print("WARNING: Cannot achieve target deHaller's ratio at hub - w2_hub would be less than cz")
    c02_hub = u1_hub  # Fallback
else:
    c02_hub = u1_hub - math.sqrt(w2_hub_sq - cz*cz)

# Calculate hub angles from c02_hub
a2_hub = math.atan(c02_hub / cz)
b2_hub = math.atan(math.tan(a2_hub) - u1_hub/cz)
b1_hub = -math.atan((u1_hub - c01) / cz)

# Verify hub deHaller's ratio
w2_hub_calc = math.sqrt(cz*cz + (u1_hub-c02_hub)*(u1_hub-c02_hub))
r1_dehallers_hub = w2_hub_calc / w1_hub

# Calculate hub enthalpy change: Δh = U * (c02 - c01)
s1dh_hub = u1_hub * (c02_hub - c01)

# ===== STEP 3: Calculate TIP angles to match midline deHaller's ratio =====
# Same process for tip

# Tip inlet relative velocity (known from geometry)
w1_tip = math.sqrt(cz*cz + (u1_tip-c01)*(u1_tip-c01))

# For tip, we need w2_tip = r1_dehallers_midline * w1_tip
w2_tip_target = r1_dehallers_midline * w1_tip

# Solve for c02_tip
w2_tip_sq = w2_tip_target * w2_tip_target
if w2_tip_sq < cz*cz:
    print("WARNING: Cannot achieve target deHaller's ratio at tip - w2_tip would be less than cz")
    c02_tip = u1_tip  # Fallback
else:
    c02_tip = u1_tip - math.sqrt(w2_tip_sq - cz*cz)

# Calculate tip angles from c02_tip
a2_tip = math.atan(c02_tip / cz)
b2_tip = math.atan(math.tan(a2_tip) - u1_tip/cz)
b1_tip = -math.atan((u1_tip - c01) / cz)

# Verify tip deHaller's ratio
w2_tip_calc = math.sqrt(cz*cz + (u1_tip-c02_tip)*(u1_tip-c02_tip))
r1_dehallers_tip = w2_tip_calc / w1_tip

# Calculate tip enthalpy change: Δh = U * (c02 - c01)
s1dh_tip = u1_tip * (c02_tip - c01)

# ===== STEP 4: Calculate degree of reaction and other parameters =====
deg_reaction_1 = 1 - (cz/(2*u1))*(math.tan(a2)+math.tan(a1))
deg_reaction_1_hub = 1 - (cz/(2*u1_hub))*(math.tan(a2_hub)+math.tan(a1))
deg_reaction_1_tip = 1 - (cz/(2*u1_tip))*(math.tan(a2_tip)+math.tan(a1))

# ===== OUTPUT RESULTS =====
print("\n" + "="*60)
print("AXIAL COMPRESSOR STAGE - UNIFORM deHALLER's RATIO DESIGN")
print("="*60)

print("\nOperating Conditions:")
print(f"  RPM: {RPM:.0f}")
print(f"  Mass Flow: {mass_flow:.2f} kg/s")
print(f"  Axial Velocity (cz): {cz:.2f} m/s")
print(f"  Rotor Work: {Wr1:.0f} J")
print(f"  Stage Loading: {s1_dh_stageloading:.4f}")

print("\nGeometry:")
print(f"  Meanline Radius: {s1_r:.2f} mm")
print(f"  Hub Radius: {inner_rad:.2f} mm")
print(f"  Tip Radius: {outer_rad:.2f} mm")

print("\nTangential Velocities (U):")
print(f"  Midline U: {u1:.2f} m/s")
print(f"  Hub U: {u1_hub:.2f} m/s")
print(f"  Tip U: {u1_tip:.2f} m/s")

print("\n" + "-"*60)
print("ANGLES (degrees):")
print("-"*60)
print(f"  IGV angle (a1): {alpha1:.2f}")
print(f"\n  Rotor Inlet Relative Angle (b1):")
print(f"    Midline: {b1*180/math.pi:.2f}")
print(f"    Hub:     {b1_hub*180/math.pi:.2f}")
print(f"    Tip:     {b1_tip*180/math.pi:.2f}")

print(f"\n  Rotor Exit Relative Angle (b2):")
print(f"    Midline: {b2*180/math.pi:.2f}")
print(f"    Hub:     {b2_hub*180/math.pi:.2f}")
print(f"    Tip:     {b2_tip*180/math.pi:.2f}")

print(f"\n  Absolute Exit Angle (a2):")
print(f"    Midline: {a2*180/math.pi:.2f}")
print(f"    Hub:     {a2_hub*180/math.pi:.2f}")
print(f"    Tip:     {a2_tip*180/math.pi:.2f}")

print(f"\n  Stator Exit Angle (a3): {a3*180/math.pi:.2f}")

print("\n" + "-"*60)
print("deHALLER's RATIOS (W2/W1):")
print("-"*60)
print(f"  Midline: {r1_dehallers_midline:.6f}")
print(f"  Hub:     {r1_dehallers_hub:.6f}")
print(f"  Tip:     {r1_dehallers_tip:.6f}")
print(f"\n  Difference from midline:")
print(f"    Hub: {abs(r1_dehallers_hub - r1_dehallers_midline):.6e}")
print(f"    Tip: {abs(r1_dehallers_tip - r1_dehallers_midline):.6e}")

print("\n" + "-"*60)
print("ENTHALPY CHANGE (Δh = U * (c02 - c01)):")
print("-"*60)
print(f"  Midline: {s1dh_midline:.2f} J/kg")
print(f"  Hub:     {s1dh_hub:.2f} J/kg")
print(f"  Tip:     {s1dh_tip:.2f} J/kg")
print(f"\n  Note: Input s1dh = {s1dh:.2f} J/kg (used for initial calculation)")

print("\n" + "-"*60)
print("DEGREE OF REACTION:")
print("-"*60)
print(f"  Midline: {deg_reaction_1:.4f}")
print(f"  Hub:     {deg_reaction_1_hub:.4f}")
print(f"  Tip:     {deg_reaction_1_tip:.4f}")

print("\n" + "-"*60)
print("TANGENTIAL VELOCITIES:")
print("-"*60)
print(f"  c01 (inlet): {c01:.2f} m/s")
print(f"  c02 (exit):")
print(f"    Midline: {c02:.2f} m/s")
print(f"    Hub:     {c02_hub:.2f} m/s")
print(f"    Tip:     {c02_tip:.2f} m/s")

print("\n" + "-"*60)
print("RELATIVE VELOCITIES:")
print("-"*60)
print(f"  W1 (inlet):")
print(f"    Midline: {w1_midline:.2f} m/s")
print(f"    Hub:     {w1_hub:.2f} m/s")
print(f"    Tip:     {w1_tip:.2f} m/s")
print(f"\n  W2 (exit):")
print(f"    Midline: {w2_midline:.2f} m/s")
print(f"    Hub:     {w2_hub_calc:.2f} m/s")
print(f"    Tip:     {w2_tip_calc:.2f} m/s")

print("\n" + "="*60)

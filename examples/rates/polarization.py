import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
import matplotlib.animation as animation

# --- Configuration ---
NUM_PARTICLES = 32
RADIUS = 2.0
AMPLITUDE = 0.8  # Strength of the GW
OMEGA = 2.0      # Frequency of oscillation

# --- Setup Particles ---
# Create a ring of particles in the Source/Sky Frame
theta = np.linspace(0, 2*np.pi, NUM_PARTICLES, endpoint=False)
x0 = RADIUS * np.cos(theta)
y0 = RADIUS * np.sin(theta)
particles_initial = np.column_stack((x0, y0))

# --- Figure Setup ---
fig, ax = plt.subplots(figsize=(8, 8))
plt.subplots_adjust(bottom=0.25, left=0.25) # Make room for controls

ax.set_xlim(-4, 4)
ax.set_ylim(-4, 4)
ax.set_aspect('equal')
ax.set_title("GW Polarization in Sky Plane\n(Propagation out of screen)")
ax.set_xlabel("Observer East (X)")
ax.set_ylabel("Observer North (Y)")
ax.grid(True, linestyle=':', alpha=0.6)

# Fixed Observer Axes (Gray)
ax.axhline(0, color='gray', linestyle='--', linewidth=0.8, label="Observer Horizon")
ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)

# --- Plot Elements ---
# 1. Particles
scatter, = ax.plot([], [], 'o', color='dodgerblue', markersize=6)

# 2. Wave Basis Vectors (The rotating axes)
# Wave X-axis (Red)
arrow_wx = ax.quiver(0, 0, 1, 0, color='red', scale=1, scale_units='xy', angles='xy', width=0.015, label="Wave X-axis")
# Wave Y-axis (Green)
arrow_wy = ax.quiver(0, 0, 0, 1, color='green', scale=1, scale_units='xy', angles='xy', width=0.015, label="Wave Y-axis")

# Text annotations for axes
text_wx = ax.text(0, 0, "Wx", color='red', fontweight='bold')
text_wy = ax.text(0, 0, "Wy", color='green', fontweight='bold')

ax.legend(loc='upper right')

# --- Interactive Controls ---

# Slider for Psi (Polarization Angle)
ax_psi = plt.axes([0.25, 0.1, 0.5, 0.03], facecolor='lightgoldenrodyellow')
slider_psi = Slider(ax_psi, r'Polarization $\psi$ (rad)', 0.0, np.pi, valinit=0.0)

# Radio Buttons for Mode Selection
ax_radio = plt.axes([0.02, 0.4, 0.15, 0.15], facecolor='#f0f0f0')
radio = RadioButtons(ax_radio, ('Plus (+)', 'Cross (x)', 'Circular'))

# --- Physics Engine ---

def get_strain(t, mode):
    """Calculates h+ and hx based on selected mode."""
    phase = OMEGA * t
    if mode == 'Plus (+)':
        hp = AMPLITUDE * np.sin(phase)
        hx = 0
    elif mode == 'Cross (x)':
        hp = 0
        hx = AMPLITUDE * np.sin(phase)
    else: # Circular
        hp = AMPLITUDE * np.cos(phase)
        hx = AMPLITUDE * np.sin(phase)
    return hp, hx

def update(frame):
    # 1. Get current params
    t = frame * 0.05
    psi = slider_psi.val
    mode = radio.value_selected
    
    # 2. Calculate Wave Strain in WAVE FRAME
    hp, hx = get_strain(t, mode)
    
    # 3. Rotate Particles: Observer Frame -> Wave Frame
    # Rotation Matrix R(-psi)
    c, s = np.cos(psi), np.sin(psi)
    
    # Project positions into Wave Frame
    x_wave =  c * x0 + s * y0
    y_wave = -s * x0 + c * y0
    
    # 4. Apply GW Distortion (Linearized Gravity)
    # dx = 0.5 * (h+ * x + hx * y)
    # dy = 0.5 * (hx * x - h+ * y)
    dx = 0.5 * (hp * x_wave + hx * y_wave)
    dy = 0.5 * (hx * x_wave - hp * y_wave)
    
    x_wave_deformed = x_wave + dx
    y_wave_deformed = y_wave + dy
    
    # 5. Rotate Back: Wave Frame -> Observer Frame
    # Rotation Matrix R(+psi)
    x_final = c * x_wave_deformed - s * y_wave_deformed
    y_final = s * x_wave_deformed + c * y_wave_deformed
    
    # 6. Update Plot
    scatter.set_data(x_final, y_final)
    
    # Update Wave Basis Vectors (Show the rotation psi)
    # Length of arrows
    L = 2.5
    arrow_wx.set_UVC(L*c, L*s)
    arrow_wy.set_UVC(-L*s, L*c)
    
    text_wx.set_position((L*c*1.1, L*s*1.1))
    text_wy.set_position((-L*s*1.1, L*c*1.1))
    
    return scatter, arrow_wx, arrow_wy, text_wx, text_wy

# --- Run Animation ---
ani = animation.FuncAnimation(fig, update, frames=200, interval=50, blit=True)

plt.show()
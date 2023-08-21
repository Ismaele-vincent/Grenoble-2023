import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.special import jv
from matplotlib.gridspec import GridSpec
matplotlib.rcParams['font.size'] = 13

mu_N=-9.6623651#*1e-27 J/T
hbar= 6.62607015/(2*np.pi) #*1e-34 J s
v0=2060.43 #m/s
phi1=0
phi2=0
order=4

# The parametrized function to be plotted
def alpha(T,f,B):
    w=f*2*np.pi
    return mu_N*B/(hbar*w)*2*np.sin(w*T*1e-3/2)

t = np.linspace(0, 1000, 1000)

# Define initial parameters
T_0=19.4 #\mu s
f1_0=1 #kHz
B1_0=10 #Gauss
f2_0=1
B2_0=0
xi1_0=phi1+(2*np.pi*f1_0*1e-3*T_0+np.pi)/2-2*np.pi*f1_0*1e3/v0
xi2_0=phi2+(2*np.pi*f2_0*1e-3*T_0+np.pi)/2-2*np.pi*f2_0*1e3/v0
chi_0=np.pi/4

a1=alpha(T_0,f1_0, B1_0)
a2=alpha(T_0,f2_0, B2_0)

# Create the figure and the line that we will manipulate
fig = plt.figure(figsize=(16,9))
title=fig.suptitle("$\\alpha1=$"+str("%.3f"%(alpha(T_0, f1_0, B1_0),))+"\n$\\alpha2=$"+str("%.3f"%(alpha(T_0, f2_0, B2_0),)))
gst = GridSpec(2,3, figure=fig, hspace=0, wspace=0)#,hspace=0, bottom=0,top=0)
gsb = GridSpec(2,3, figure=fig, hspace=0, wspace=0)
axs = [fig.add_subplot(gst[:,:-1]),
       fig.add_subplot(gst[0,-1]),
       fig.add_subplot(gsb[1,-1])]

for ax in axs[1:]:
    ax.yaxis.set_tick_params(labelleft=False)
    ax.yaxis.set_label_position("right")
    
axs[0].set_xlabel("Time ($\mu$s)")
axs[0].set_ylabel("$I_o$")

# adjust the main plot to make room for the sliders
lineO, = axs[0].plot(t,(1+np.cos(chi_0-a1*np.sin(2*np.pi*1e-3*f1_0*t+xi1_0)+a2*np.sin(2*np.pi*1e-3*f2_0*t+xi2_0)))/2, "k-", ms=2)

axs[1].set_xlabel("Order")
axs[1].set_ylabel("$J_n^2$ Path1")
rects1=axs[1].bar(range(-order, order+1), jv(range(-order, order+1), alpha(T_0, f1_0, B1_0))**2)
texts1=[]
for i in range(-order,order+1):
    text=axs[1].text(i, jv(i, alpha(T_0, f1_0, B1_0))**2, str("%.3f" %(jv(i, alpha(T_0, f1_0, B1_0))**2),),ha="center", va="bottom")
    texts1.append(text)
axs[1].set_xticks(range(-order, order+1))
axs[1].set_ylim([0,1.1])    

axs[2].set_xlabel("Order")
axs[2].set_ylabel("$J_n^2$ Path2")
rects2=axs[2].bar(range(-order, order+1), jv(range(-order, order+1), alpha(T_0, f2_0, B2_0))**2)
texts2=[]
for i in range(-order,order+1):
    text=axs[2].text(i, jv(i, alpha(T_0, f2_0, B2_0))**2, str("%.3f" %(jv(i, alpha(T_0, f2_0, B2_0))**2),),ha="center", va="bottom")
    texts2.append(text)
axs[2].set_xticks(range(-order, order+1))
axs[2].set_ylim([0,1.1])

# Make horizontal sliders.
fig.subplots_adjust(bottom=0.25, left=0.08)
axf1 = fig.add_axes([0.25, 0.13, 0.5, 0.02])
axf2 = fig.add_axes([0.25, 0.10, 0.5, 0.02])
axB1 = fig.add_axes([0.25, 0.07, 0.5, 0.02])
axB2 = fig.add_axes([0.25, 0.04, 0.5, 0.02])
axT = fig.add_axes([0.25, 0.01, 0.5, 0.02])
axchi = fig.add_axes([0.01, 0.25, 0.02, 0.5])

f1_slider = Slider(
    ax=axf1,
    label="$f_1$ (kHz)",
    valmin=0.1,
    valmax=80,
    valinit=f1_0,
)

f2_slider = Slider(
    ax=axf2,
    label="$f_2$ (kHz)",
    valmin=0.1,
    valmax=80,
    valinit=f2_0,
)


B1_slider = Slider(
    ax=axB1,
    label="$B_1$ (Gauss)",
    valmin=0,
    valmax=60,
    valinit=B1_0,
    # orientation="vertical"
)


B2_slider = Slider(
    ax=axB2,
    label="$B_2$ (Gauss)",
    valmin=0,
    valmax=60,
    valinit=B2_0,
    # orientation="vertical"
)

T_slider = Slider(
    ax=axT,
    label="T ($\mu$s)",
    valmin=0,
    valmax=50,
    valinit=T_0,
    # orientation="vertical"
)

chi_slider = Slider(
    ax=axchi,
    label="$\chi$",
    valmin=0,
    valmax=2*np.pi,
    valinit=chi_0,
    orientation="vertical"
)

# The function to be called anytime a slider"s value changes
def update(val):
    title.set_text("$\\alpha=$"+str("%.3f"%(alpha(T_slider.val, f1_slider.val, B1_slider.val),))+"\n$\\alpha2=$"+str("%.3f"%(alpha(T_slider.val, f2_slider.val, B2_slider.val),)))
    xi1=phi1+(2*np.pi*f1_slider.val*1e-3*T_slider.val+np.pi)/2-2*np.pi*f1_slider.val*1e3/(v0)
    xi2=phi2+(2*np.pi*f2_slider.val*1e-3*T_slider.val+np.pi)/2-2*np.pi*f2_slider.val*1e3/(v0)
    a1=alpha(T_slider.val, f1_slider.val, B1_slider.val)
    a2=alpha(T_slider.val, f2_slider.val, B2_slider.val)
    lineO.set_ydata((1+np.cos(chi_slider.val-a1*np.sin(2*np.pi*f1_slider.val*1e-3*t+xi1)+a2*np.sin(2*np.pi*1e-3*f2_slider.val*t+xi2)))/2)
    for rect, h in zip(rects1, jv(range(-order, order+1), alpha(T_slider.val, f1_slider.val, B1_slider.val))**2):
        rect.set_height(h)
    for text,i in zip(texts1,range(-order,order+1)):
        text.set_text(str("%.3f" %(jv(i, alpha(T_slider.val, f1_slider.val, B1_slider.val))**2),))
        text.set_y((jv(i, alpha(T_slider.val, f1_slider.val, B1_slider.val))**2))
    for rect, h in zip(rects2, jv(range(-order, order+1), alpha(T_slider.val, f2_slider.val, B2_slider.val))**2):
        rect.set_height(h)
    for text,i in zip(texts2,range(-order,order+1)):
        text.set_text(str("%.3f" %(jv(i, alpha(T_slider.val, f2_slider.val, B2_slider.val))**2),))
        text.set_y((jv(i, alpha(T_slider.val, f2_slider.val, B2_slider.val))**2))    
        axs[0].relim()
        axs[0].autoscale_view()
    fig.canvas.draw_idle()

# register the update function with each slider
f1_slider.on_changed(update)
f2_slider.on_changed(update)
B1_slider.on_changed(update)
B2_slider.on_changed(update)
chi_slider.on_changed(update)
T_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = fig.add_axes([0.85, 0.025, 0.1, 0.04])
buttonr = Button(resetax, "Reset", hovercolor="0.975")
def reset(event):
    f1_slider.reset()
    f2_slider.reset()
    B1_slider.reset()
    B2_slider.reset()
    T_slider.reset()
    chi_slider.reset()
buttonr.on_clicked(reset)

setxi = fig.add_axes([0.85, 0.1, 0.1, 0.04])
buttons = Button(setxi, "Set $\\xi_i=\pi$", hovercolor="0.975")
def setxi(event):
    xi1=np.pi
    xi2=np.pi
    a1=alpha(T_slider.val, f1_slider.val, B1_slider.val)
    a2=alpha(T_slider.val, f2_slider.val, B2_slider.val)
    lineO.set_ydata((1+np.cos(chi_slider.val-a1*np.sin(2*np.pi*f1_slider.val*1e-3*t+xi1)+a2*np.sin(2*np.pi*1e-3*f2_slider.val*t+xi2)))/2)
buttons.on_clicked(setxi)

plt.show()
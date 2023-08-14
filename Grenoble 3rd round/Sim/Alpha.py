import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.special import jv
from matplotlib.gridspec import GridSpec

mu_N=-9.6623651#*1e-27 J/T
hbar= 6.62607015/(2*np.pi) #*1e-34 J s
order=3

# The parametrized function to be plotted
def alpha(T,f1,B1):
    w1=f1*2*np.pi
    B1*=1e-4
    return mu_N*B1*1e4/(hbar*w1)*2*np.sin(w1*T*1e-3/2)



T0=19.4 #\mu s
T = np.linspace(0, 100, 10000)

# Define initial parameters
f1_0=1 #kHz
B10=10 #Gauss

# Create the figure and the line that we will manipulate
fig = plt.figure(figsize=(16,9))
gst = GridSpec(2,2, figure=fig, hspace=0, bottom=0.2)#,hspace=0, bottom=0,top=0)
gsb = GridSpec(2,2, figure=fig, hspace=0, bottom=0.2)
axs = [fig.add_subplot(gst[0,0]),
       fig.add_subplot(gst[1,0]),
       fig.add_subplot(gsb[:,1])]
line, = axs[0].plot(T, alpha(T, f1_0, B10), lw=2)
for ax in axs:
    ax.yaxis.set_label_position("right")
axs[0].set_xlabel("Time ($\mu$ s)")
axs[0].set_ylabel("$\\alpha$")
fig.suptitle("$T_0$="+str(T0)+" $\mu s$\n$\\alpha=$"+str("%.3f"%(alpha(T0, f1_0, B10),)))
axs[0].axvline(T0, ls="dashed", color="k", label="$T_0$")
axs[0].legend(loc=1, framealpha=1)
axs[1].set_xlabel("Time ($\mu$ s)")
axs[1].set_ylabel("$J_n(\\alpha)^2$")
axs[1].axvline(T0, ls="dashed", color="k")
# adjust the main plot to make room for the sliders
lineJ0, = axs[1].plot(T, jv(0,alpha(T, f1_0, B10))**2, "k-", ms=2, label="$J_0^2$")
lineJ1p, = axs[1].plot(T, jv(1,alpha(T, f1_0, B10))**2, "-", ms=2, label="$J_1^2$")
lineJ2p, = axs[1].plot(T, jv(2,alpha(T, f1_0, B10))**2, "-", ms=2, label="$J_{2}^2$")
lineJ3p, = axs[1].plot(T, jv(3,alpha(T, f1_0, B10))**2, "-", ms=2, label="$J_{3}^2$")
lineJ4p, = axs[1].plot(T, jv(4,alpha(T, f1_0, B10))**2, "-", ms=2, label="$J_{4}^2$")

axs[1].legend(loc=1)
axs[2].set_xlabel("Order")
axs[2].set_ylabel("$J_n^2 (T=T_0)$")
axs[2].bar(range(-order, order+1), jv(range(-order, order+1), alpha(T0, f1_0, B10))**2)
axs[2].text(0, jv(0, alpha(T0, f1_0, B10))**2, str("%.3f" %(jv(0, alpha(T0, f1_0, B10))**2),),ha="center", va="bottom")
for i in range(1,order+1):
    axs[2].text(i, jv(i, alpha(T0, f1_0, B10))**2, str("%.3f" %(jv(i, alpha(T0, f1_0, B10))**2),),ha="center", va="bottom")
    axs[2].text(-i, jv(i, alpha(T0, f1_0, B10))**2, str("%.3f" %(jv(-i, alpha(T0, f1_0, B10))**2),),ha="center", va="bottom")
axs[2].set_xticks(range(-order, order+1))
axs[2].set_ylim([0,1])
# Make horizontal sliders.
axB1 = fig.add_axes([0.25, 0.05, 0.5, 0.03])
fig.subplots_adjust(bottom=0.1)
axf1 = fig.add_axes([0.25, 0.01, 0.5, 0.03])
f1_slider = Slider(
    ax=axf1,
    label="$f_1$ (kHz)",
    valmin=0.1,
    valmax=80,
    valinit=f1_0,
)

B1_slider = Slider(
    ax=axB1,
    label="B1 (Gauss)",
    valmin=0,
    valmax=60,
    valinit=B10,
    # orientation="vertical"
)


# The function to be called anytime a slider"s value changes
def update(val):
    line.set_ydata(alpha(T, f1_slider.val, B1_slider.val))
    fig.suptitle("$T_0$="+str(T0)+" $\mu s$\n$\\alpha=$"+str("%.3f"%(alpha(T0, f1_slider.val, B1_slider.val),)))
    lineJ0.set_ydata(jv(0,alpha(T, f1_slider.val, B1_slider.val))**2)
    lineJ1p.set_ydata(jv(1,alpha(T, f1_slider.val, B1_slider.val))**2)
    lineJ2p.set_ydata(jv(2,alpha(T, f1_slider.val, B1_slider.val))**2)
    lineJ3p.set_ydata(jv(3,alpha(T, f1_slider.val, B1_slider.val))**2)
    lineJ4p.set_ydata(jv(4,alpha(T, f1_slider.val, B1_slider.val))**2)
    axs[2].clear()
    axs[2].bar(range(-order, order+1), jv(range(-order, order+1), alpha(T0, f1_slider.val, B1_slider.val))**2)
    axs[2].set_xlabel("Order")
    axs[2].set_xticks(range(-order, order+1))
    axs[2].set_ylim([0,1])
    axs[2].text(0, jv(0, alpha(T0, f1_slider.val, B1_slider.val))**2, str("%.3f" %(jv(0, alpha(T0, f1_slider.val, B1_slider.val))**2),),ha="center", va="bottom")
    for i in range(1,order+1):
        axs[2].text(i, jv(i, alpha(T0, f1_slider.val, B1_slider.val))**2, str("%.3f" %(jv(i, alpha(T0, f1_slider.val, B1_slider.val))**2),),ha="center", va="bottom")
        axs[2].text(-i, jv(i, alpha(T0, f1_slider.val, B1_slider.val))**2, str("%.3f" %(jv(i, alpha(T0, f1_slider.val, B1_slider.val))**2),),ha="center", va="bottom")
    for ax in axs[:-2]:    
        ax.relim()
        ax.autoscale_view()
    fig.canvas.draw_idle()


# register the update function with each slider
f1_slider.on_changed(update)
B1_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = fig.add_axes([0.85, 0.025, 0.1, 0.04])
button = Button(resetax, "Reset", hovercolor="0.975")


def reset(event):
    f1_slider.reset()
    B1_slider.reset()
button.on_clicked(reset)

plt.show()






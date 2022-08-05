from grlp import *
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap

# ------------------------------------------------------------------------------

# Figure x in McNab et al.

# ------------------------------------------------------------------------------


# ---- River properties
L = 100.e3
Qw = 16.
B = 150.
dx = 1.e3
Qs0 = 0.001163080432007373

# ---- Set up long profile object
lp = LongProfile()
lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()
lp.set_x(x_ext=np.arange(0., L+dx, dx))
lp.set_Q(Q=Qw)
lp.set_B(B=B)
lp.set_z(S0=(Qs0/(lp.k_Qs*Qw))**(6./7.))
lp.set_Qs_input_upstream(Qs0)
lp.set_z_bl(0.)
lp.set_uplift_rate(0.)
lp.set_niter()
lp.evolve_threshold_width_river(nt=1000, dt=3.15e10)
lp.compute_Q_s()
lp.compute_equilibration_time()

# ---- Record properties
lin_periods = np.logspace(-2., 2., 81) * 3.15e12
lin_gain = np.zeros((len(lin_periods), len(lp.x)))
lin_lag = np.zeros((len(lin_periods), len(lp.x)))
lin_gain_Qs = np.zeros((len(lin_periods), len(lp.x)))
lin_lag_Qs = np.zeros((len(lin_periods), len(lp.x)))
lin_gain_Qs_Qw = np.zeros((len(lin_periods), len(lp.x)))
lin_lag_Qs_Qw = np.zeros((len(lin_periods), len(lp.x)))
for i,p in enumerate(lin_periods):
    lin_gain[i,:] = lp.compute_z_gain(p)
    lin_lag[i,:] = lp.compute_z_lag(p) / p
    lin_gain_Qs[i,:] = lp.compute_Qs_gain(p, A_Qs=0.2)
    lin_lag_Qs[i,:] = lp.compute_Qs_lag(p, A_Qs=0.2) / p
    lin_gain_Qs_Qw[i,:] = lp.compute_Qs_gain(p, A_Q=0.2)
    lin_lag_Qs_Qw[i,:] = lp.compute_Qs_lag(p, A_Q=0.2) / p

# ---- L'
L_dash = np.sqrt(lin_periods * lp.diffusivity.mean())

# ---- Plot grids

# get colours & contours for gain
red_to_black = LinearSegmentedColormap.from_list("", ["tomato", "black"])
blue_to_red = LinearSegmentedColormap.from_list("", ["cornflowerblue", "white", "tomato"])
gain_color_lvls = np.arange(0., 1.201, 0.001)
gain_contour_lvls = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2]
gain_colors = np.vstack(( 
    blue_to_red(gain_color_lvls[gain_color_lvls < 1]), 
    red_to_black((gain_color_lvls[gain_color_lvls >= 1] - 1.) / \
        max(gain_color_lvls[gain_color_lvls >= 1] - 1.))
    ))
gain_cbar_lvls = [0., 0.2, 0.4, 0.6, 0.8, 1., 1.2]

# get colours & contours for lag
hot = cm.get_cmap("hot_r", 100)
bone = cm.get_cmap("bone_r", 100)
lag_color_lvls1 = np.arange(0., 3.01, 0.01)
lag_contour_lvls1 = np.arange(0., 3., 0.25)
lag_colors1 = hot((lag_color_lvls1)/(lag_color_lvls1+0.2).max())
lag_color_lvls2 = np.arange(-0.25, 0., 0.01)
lag_contour_lvls2 = np.arange(-0.25, 0., 0.05)
lag_colors2 = bone((lag_color_lvls2+0.07)/lag_color_lvls2.min())
lag_color_lvls = np.hstack(( lag_color_lvls2, lag_color_lvls1 ))
lag_contour_lvls = np.hstack(( lag_contour_lvls2, lag_contour_lvls1))
lag_colors = np.vstack(( lag_colors2, lag_colors1, ))
lag_cbar_lvls = [-0.25, 0., 0.5, 1., 1.5, 2., 2.5, 3.]

# initialise plot
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(12,7))

# dz, gain
ax1.contourf(lin_periods/3.15e12, lp.x/lp.L, lin_gain.swapaxes(0,1), levels=gain_color_lvls, colors=gain_colors)
ax1.contour(lin_periods/3.15e12, lp.x/lp.L, lin_gain.swapaxes(0,1), colors="k", linewidths=0.5, levels=gain_contour_lvls)
ax1.plot(lin_periods/3.15e12, L_dash/L, "--", linewidth=1.5, color="k")
label = r"$x = \sqrt{P\kappa}$"
ax1.annotate(label, (lin_periods[20]/3.15e12, L_dash[20]/L), ha="right")
ax1.set_ylabel(r"$x~/~L$")
ax1.set_box_aspect(1)
ax1.set_ylim(0,1)
ax1.set_title(r"${\delta} z$")

# dQ_s, varying sediment supply, gain
ax2.contourf(lin_periods/3.15e12, lp.x/lp.L, lin_gain_Qs.swapaxes(0,1), levels=gain_color_lvls, colors=gain_colors)
ax2.contour(lin_periods/3.15e12, lp.x/lp.L, lin_gain_Qs.swapaxes(0,1), colors="k", linewidths=0.5, levels=gain_contour_lvls)
ax2.plot(lin_periods/3.15e12, L_dash/L, "--", linewidth=1.5, color="k")
ax2.set_box_aspect(1)
ax2.set_title(r"${\delta} Q_s$: Varying sediment supply")

# dQ_s, varying water supply, gain
cmap = ax3.contourf(lin_periods/3.15e12, lp.x/lp.L, lin_gain_Qs_Qw.swapaxes(0,1), levels=gain_color_lvls, colors=gain_colors)
ax3.contour(lin_periods/3.15e12, lp.x/lp.L, lin_gain_Qs_Qw.swapaxes(0,1), colors="k", linewidths=0.5, levels=gain_contour_lvls)
ax3.plot(lin_periods/3.15e12, L_dash/L, "--", linewidth=1.5, color="k")
ax3.set_box_aspect(1)
cax = ax3.inset_axes([1.1, 0., 0.05, 1.], transform=ax3.transAxes)
cbar = fig.colorbar(cmap, ax=ax3, cax=cax)
cbar.set_ticks(gain_cbar_lvls)
cbar.set_ticklabels(gain_cbar_lvls)
cbar.set_label(r"Gain, $G$")
ax3.set_title(r"${\delta} Q_s$: Varying water supply")

# dz, lag
ax4.contourf(lin_periods/3.15e12, lp.x/lp.L, lin_lag.swapaxes(0,1), levels=lag_color_lvls, colors=lag_colors)
ax4.contour(lin_periods/3.15e12, lp.x/lp.L, lin_lag.swapaxes(0,1), colors="k", linewidths=0.5, levels=lag_contour_lvls)
ax4.plot(lin_periods/3.15e12, L_dash/L, "--", linewidth=1.5, color="k")
ax4.set_xscale("log")
ax4.set_xlabel(r"$P~/~T_{eq}$")
ax4.set_ylabel(r"$x~/~L$")
ax4.set_box_aspect(1)

# dQ_s, varying sediment supply, lag
ax5.contourf(lin_periods/3.15e12, lp.x/lp.L, lin_lag_Qs.swapaxes(0,1), levels=lag_color_lvls, colors=lag_colors)
ax5.contour(lin_periods/3.15e12, lp.x/lp.L, lin_lag_Qs.swapaxes(0,1), colors="k", linewidths=0.5, levels=lag_contour_lvls)
ax5.plot(lin_periods/3.15e12, L_dash/L, "--", linewidth=1.5, color="k")
ax5.set_xlabel(r"$P~/~T_{eq}$")
ax5.set_box_aspect(1)

# dQ_s, varying water supply, lag
cmap_l = ax6.contourf(lin_periods/3.15e12, lp.x/lp.L, lin_lag_Qs_Qw.swapaxes(0,1), levels=lag_color_lvls, colors=lag_colors)
ax6.contour(lin_periods/3.15e12, lp.x/lp.L, lin_lag_Qs_Qw.swapaxes(0,1), colors="k", linewidths=0.5, levels=lag_contour_lvls)
ax6.plot(lin_periods/3.15e12, L_dash/L, "--", linewidth=1.5, color="k")
ax6.set_xlabel(r"$P~/~T_{eq}$")
ax6.set_box_aspect(1)
cax = ax6.inset_axes([1.1, 0., 0.05, 1.], transform=ax6.transAxes)
cbar = fig.colorbar(cmap_l, ax=ax6, cax=cax)
cbar.set_ticks(lag_cbar_lvls)
cbar.set_ticklabels(lag_cbar_lvls)
cbar.set_label(r"Lag, $\varphi~/~P$")

plt.tight_layout()
plt.show()



# # ---- Save data for GMT
# # Various files used to make a nice version of the figure in GMT. 
#
# z_out = "../../../periodic_output_1D/sweeps/"
# Qs_out = "../../../periodic_output_1D/Qs_sweeps/"
# 
# with open(z_out + "linear.pdg", "wb") as f:
#     for i,p in enumerate(lin_periods):
#         arr = np.column_stack((
#             np.full(len(lp.x), p/lp.equilibration_time),
#             lp.x/1000.,
#             lin_gain[i,:]
#             ))
#         np.savetxt(f, arr)
# 
# with open(z_out + "linear.pdl", "wb") as f:
#     for i,p in enumerate(lin_periods):
#         arr = np.column_stack((
#             np.full(len(lp.x), p/lp.equilibration_time),
#             lp.x/1000.,
#             lin_lag[i,:]
#             ))
#         np.savetxt(f, arr)
# 
# with open(Qs_out + "linear_Qs.pdg", "wb") as f:
#     for i,p in enumerate(lin_periods):
#         arr = np.column_stack((
#             np.full(len(lp.x), p/lp.equilibration_time),
#             lp.x/1000.,
#             lin_gain_Qs[i,:]
#             ))
#         np.savetxt(f, arr)
# 
# with open(Qs_out + "linear_Qs.pdl", "wb") as f:
#     for i,p in enumerate(lin_periods):
#         arr = np.column_stack((
#             np.full(len(lp.x), p/lp.equilibration_time),
#             lp.x/1000.,
#             lin_lag_Qs[i,:]
#             ))
#         np.savetxt(f, arr)
# 
# 
# with open(Qs_out + "linear_Qw.pdg", "wb") as f:
#     for i,p in enumerate(lin_periods):
#         arr = np.column_stack((
#             np.full(len(lp.x), p/lp.equilibration_time),
#             lp.x/1000.,
#             lin_gain_Qs_Qw[i,:]
#             ))
#         np.savetxt(f, arr)
# 
# with open(Qs_out + "linear_Qw.pdl", "wb") as f:
#     for i,p in enumerate(lin_periods):
#         arr = np.column_stack((
#             np.full(len(lp.x), p/lp.equilibration_time),
#             lp.x/1000.,
#             lin_lag_Qs_Qw[i,:]
#             ))
#         np.savetxt(f, arr)
# 
# with open(Qs_out + "L_dash_contour.pd", "wb") as f:
#     arr = np.column_stack((
#         lin_periods/lp.equilibration_time, 
#         np.sqrt(lin_periods*lp.diffusivity.mean())/1000./100. ))
#     np.savetxt(f, arr)
from grlp import *

def predict_hack_length(area, k=1.85, h=0.54):
    return k*(area**h)

def predict_hack_area(length, k=1.85, h=0.54):
    return (length/k)**(1./h)
    
def predict_mean_hack_area(length, k=1.85, h=0.54, x0=10.e3):
    return (((length+x0)/k)**(1./h + 1.) - (x0/k)**(1./h + 1.)) / (length*(1./h + 1))
    
def predict_clubb_width(area, k, p):
    return k * (area**p)

def predict_mean_clubb_width(area, k_b, p_b):
    return k_b * (area**p_b) / (p_b + 1.)

def predict_mean_clubb_width_hack(length, k_b, p_b, k=1.85, h=0.54, x0=10.e3):
    return k_b * ((length+x0)**(p_b/h + 1.) - x0**(p_b/h + 1.)) / ((length+x0) * (k**(p_b/h)) * (p_b/h + 1.))
    
def predict_precipitation(length, T_eq, k_b=0.4687028725672903, p_b=0.3, slope=0.05, k_Qs=0.040987384904837776, x0=10.e3, runoff_coeff=0.4):
    meanA = predict_mean_hack_area(length, x0=x0)
    # meanL = predict_hack_length(meanA)
    meanB = predict_mean_clubb_width_hack(length, k_b, p_b, x0=x0)
    # meanB = predict_mean_clubb_width(predict_hack_area(length), k_b, p_b, x0=x0)
    # meanB = 100.
    return (length)**2. * meanB * (1. - 0.35) / ((7./6.) * T_eq * runoff_coeff * k_Qs * meanA * slope**(1./6.))


lengths = np.linspace(1.e3, 200.e3, 100)
periods = np.array([20.,40.,100.,400.]) * 3.154e10
slopes = np.array([0.002, 0.01, 0.05])
P_steep = np.zeros((len(lengths),len(periods)))
P_med = np.zeros((len(lengths),len(periods)))
P_shallow = np.zeros((len(lengths),len(periods)))

for i,p in enumerate(periods):
    P_steep[:,i] = predict_precipitation(lengths, p, slope=slopes.max())
    P_med[:,i] = predict_precipitation(lengths, p, slope=slopes[1])
    P_shallow[:,i] = predict_precipitation(lengths, p, slope=slopes.min())

secs_in_year = 3.154e7
mm_in_metre = 1.e3

for i,p in enumerate(periods):
    plt.fill(
        np.hstack((lengths,lengths[::-1]))/1.e3,
        np.hstack((P_steep[:,i],P_shallow[:,i][::-1]))*secs_in_year*mm_in_metre,
        alpha=0.2
        )
    lab = r"$T_{eq}$ = %i kyr" % (p/3.15e10)
    plt.plot(lengths/1.e3, P_med[:,i]*secs_in_year*mm_in_metre, label=lab)
plt.xlim(0., 200.)
plt.ylim(1.e2, 1.e4)
plt.yscale("log")
plt.xlabel(r"Valley length, $L$ [km]")
plt.ylabel(r"Precipitation rate [mm yr$^{-1}$]")
plt.legend()
plt.gca().set_aspect(1./plt.gca().get_data_ratio())
plt.show()


with open("./figures/eq_times.dat", "wb") as f:
    for i,p in enumerate(periods):
        hdr = b'> -L"%i kyr" -Z%f\n' % (p/3.15e10, p/3.15e10)
        f.write(hdr)
        arr = np.column_stack(( lengths/1.e3, P_med[:,i]*secs_in_year*mm_in_metre ))
        np.savetxt(f, arr)

with open("./figures/eq_times_poly.dat", "wb") as f:
    for i,p in enumerate(periods):
        hdr = b"> -Z%f\n" % (p/3.15e10)
        f.write(hdr)
        arr = np.column_stack(( lengths/1.e3, P_steep[:,i]*secs_in_year*mm_in_metre ))
        np.savetxt(f, arr)
        arr = np.column_stack(( lengths[::-1]/1.e3, P_shallow[:,i][::-1]*secs_in_year*mm_in_metre ))
        np.savetxt(f, arr)


sys.exit()






lp = LongProfile()
lp.basic_constants()
lp.bedload_lumped_constants()
k_Qs = lp.k_Qs

L_min = 1.e3
L_max = 200.e3

hack_p = 0.6
hack_k = 0.6


L = np.linspace(L_min, L_max, 100)
meanA = predict_mean_hack_area(L)

T = np.array([20., 200., 40., 400., 100., 1000.])*3.15e10
P = np.zeros((len(T), len(L)))
P_min = np.zeros((len(T), len(L)))
P_max = np.zeros((len(T), len(L)))
for i,t in enumerate(T):
    P[i,:] = predict_precipitation(L, t)
    # P_min[i,:] = predict_precipitation(L, t, 0.086, 0.36)
    # P_max[i,:] = predict_precipitation(L, t, 2.072, 0.24)

# Toro
obsP_100 = np.array([5.43914e-09, 5.79211e-09])
obsA_100 = np.array([2.47872e+08, 3.62031e+08])
obsL_100 = np.array([90.e3, 152861.])

# Huerta
obsP_40 = np.array([7.39883e-09])
obsA_40 = np.array([7.5097e+07])
obsL_40 = np.array([28495.])

# Niquizanga, Agua Brava
obsP_20 = np.array([9.81137e-09, 9.3228e-09])
obsA_20 = np.array([5.35948e+07, 1.74178e+07])
obsL_20 = np.array([39665.7, 10368.3])

# # plot
# fig, (ax1, ax2, ax3) = plt.subplots(1,3,sharey=True,sharex=True)
# 
# # 20 kyr
# ax1.plot(meanA/1.e6, P[0,:]*3.154e+10)
# ax1.plot(meanA/1.e6, P_min[0,:]*3.154e+10, ":")
# ax1.plot(meanA/1.e6, P_max[0,:]*3.154e+10, ":")
# ax1.plot(meanA/1.e6, P[1,:]*3.154e+10)
# ax1.plot(meanA/1.e6, P_min[1,:]*3.154e+10, ":")
# ax1.plot(meanA/1.e6, P_max[1,:]*3.154e+10, ":")
# 
# # 40 kyr
# ax2.plot(meanA/1.e6, P[2,:]*3.154e+10)
# ax2.plot(meanA/1.e6, P_min[2,:]*3.154e+10, ":")
# ax2.plot(meanA/1.e6, P_max[2,:]*3.154e+10, ":")
# ax2.plot(meanA/1.e6, P[3,:]*3.154e+10)
# ax2.plot(meanA/1.e6, P_min[3,:]*3.154e+10, ":")
# ax2.plot(meanA/1.e6, P_max[3,:]*3.154e+10, ":")
# ax2.scatter(obsA_40/1.e6, obsP_40*3.154e+10)
# 
# # 100 kyr
# ax3.plot(meanA/1.e6, P[4,:]*3.154e+10)
# ax3.plot(meanA/1.e6, P_min[4,:]*3.154e+10, ":")
# ax3.plot(meanA/1.e6, P_max[4,:]*3.154e+10, ":")
# ax3.plot(meanA/1.e6, P[5,:]*3.154e+10)
# ax3.plot(meanA/1.e6, P_min[5,:]*3.154e+10, ":")
# ax3.plot(meanA/1.e6, P_max[5,:]*3.154e+10, ":")
# ax3.scatter(obsA_100/1.e6, obsP_100*3.154e+10)
# 
# ax3.set_ylim(10.,1.e4)
# ax3.set_yscale("log")
# plt.show()

# plot
fig, (ax1, ax2, ax3) = plt.subplots(1,3,sharey=True,sharex=True)

# 20 kyr
ax1.plot(L/1.e3, P[0,:]*3.154e+10)
ax1.plot(L/1.e3, P_min[0,:]*3.154e+10, ":")
ax1.plot(L/1.e3, P_max[0,:]*3.154e+10, ":")
ax1.plot(L/1.e3, P[1,:]*3.154e+10)
ax1.plot(L/1.e3, P_min[1,:]*3.154e+10, ":")
ax1.plot(L/1.e3, P_max[1,:]*3.154e+10, ":")
ax1.scatter(obsL_20/1.e3, obsP_20*3.154e+10)

# 40 kyr
ax2.plot(L/1.e3, P[2,:]*3.154e+10)
ax2.plot(L/1.e3, P_min[2,:]*3.154e+10, ":")
ax2.plot(L/1.e3, P_max[2,:]*3.154e+10, ":")
ax2.plot(L/1.e3, P[3,:]*3.154e+10)
ax2.plot(L/1.e3, P_min[3,:]*3.154e+10, ":")
ax2.plot(L/1.e3, P_max[3,:]*3.154e+10, ":")
ax2.scatter(obsL_40/1.e3, obsP_40*3.154e+10)

# 100 kyr
ax3.plot(L/1.e3, P[4,:]*3.154e+10)
ax3.plot(L/1.e3, P_min[4,:]*3.154e+10, ":")
ax3.plot(L/1.e3, P_max[4,:]*3.154e+10, ":")
ax3.plot(L/1.e3, P[5,:]*3.154e+10)
ax3.plot(L/1.e3, P_min[5,:]*3.154e+10, ":")
ax3.plot(L/1.e3, P_max[5,:]*3.154e+10, ":")
ax3.scatter(obsL_100/1.e3, obsP_100*3.154e+10)

ax3.set_ylim(10.,1.e4)
ax3.set_yscale("log")
plt.show()
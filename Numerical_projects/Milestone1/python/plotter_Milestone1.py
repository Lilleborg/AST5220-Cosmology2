import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def color_each_regime(ax,x_array,redshift=[False]):
    """ 
    Read x-values for equality between the different regimes and color in the
    background of the axes with appropriate colors.

    @ ax - axis object to fill

    @ x_array - array, the logarithmic scale factor array used for the data in figure

    @ redshift - array, default False. If filling axis plotting with redshift, give redshift array
    """
    # Reading x-values for equality from file. To avoid handling many input arguments this
    # is read each time which is not very efficient, but not so concerned about speed in the plotters.
    # Timed for 100 calls, which is way more than I ever expect to use it, which takes about 1.9 seconds.
    # 18.7 seconds for 1000 calls.
    x_values = np.loadtxt("../../data/x_values_rad_mat_lambda_eq.txt")
    idx_rad_mat_eq = np.argmin(np.abs(x_array-x_values[0]))
    idx_mat_lambda_eq = np.argmin(np.abs(x_array-x_values[1]))
    if len(redshift) < 2:
        ax.axvspan(x_array.min(), x_array[idx_rad_mat_eq], alpha=0.2, color='C0')
        ax.axvspan(x_array[idx_rad_mat_eq], x_array[idx_mat_lambda_eq], alpha=0.2, color='C1')
        ax.axvspan(x_array[idx_mat_lambda_eq], x_array.max(), alpha=0.2, color='C3')
    else:
        ax.axvspan(redshift.max(), redshift[idx_rad_mat_eq], alpha=0.2, color='C0')
        ax.axvspan(redshift[idx_rad_mat_eq], redshift[idx_mat_lambda_eq], alpha=0.2, color='C1')
        ax.axvspan(redshift[idx_mat_lambda_eq], redshift.min(), alpha=0.2, color='C3')

if __name__== "__main__":
    # Plot styling
    plt.style.use('bmh')
    mpl.rc('font', family='serif', size=15)
    mpl.rc('text', usetex=False)

    # Load data
    # ---------
    x_array, eta_of_x, H_of_x, dHdx_of_x, Hp_of_x, dHpdx_of_x, ddHpddx_of_x, OmegaB, OmegaCDM, \
        OmegaLambda, OmegaR, _, _ = np.loadtxt("../../data/cosmology.txt", unpack=True)

    # Data handling and converting
    # ----------------------------
    z_array = 1/np.exp(x_array) - 1
    z_non_neg = z_array[z_array > 0]
    a_array = np.exp(x_array)
    x_array_ticks = np.append(np.linspace(x_array.min(), x_array.max(), 5),0)
    z_array_ticks = np.logspace(7,-2,4)

    # Indices for radiation-matter and matter-lambda equalities:
    Omega_matter = OmegaB + OmegaCDM        # Omega matter combined
    idx_rad_mat_eq = np.argmin(np.abs(OmegaR-Omega_matter))
    idx_mat_lambda_eq = np.argmin(np.abs(OmegaLambda-Omega_matter))
    # Saving value of x in these indices to file to be used in plotting in later Milestones
    np.savetxt("../../data/x_values_rad_mat_lambda_eq.txt",[x_array[idx_rad_mat_eq],x_array[idx_mat_lambda_eq]])

    # Some units
    km = 1e3                                # [m]
    Mpc = 3.08567758e22                     # [m]
    h = 0.7                                 # unitless hubble
    H0 = 100 * h * km/Mpc                   # Hubble parameter today
    c = 299792458                           # Speed of light [m/s]
    # H-prime in presentable units [km/s/Mpc]:
    Hp_units = Hp_of_x * Mpc/km
    H_units = Hp_units / a_array            # H (not prime)  [km/s/Mpc]
    eta_units = eta_of_x / (Mpc*1e3)        # conformal time [Gpc]

    # Analytical etas for each regime:
    aeta_R = c/Hp_of_x/Mpc/km                      # analytical eta in radiation dominated regime
    idx_astar = np.argmax(Omega_matter)
    aeta_M = (eta_of_x[idx_astar] + 2*c*(1/Hp_of_x-1/Hp_of_x[idx_astar]))/Mpc/km
    idx_alamda = np.argmax(OmegaLambda)
    aeta_Lambda = (eta_of_x[idx_alamda] - 3*c*(1/(a_array**2*Hp_of_x) -
                                            1/(a_array[idx_alamda]**2*Hp_of_x[idx_alamda])))/Mpc/km

    # Plotting
    # --------
    # Omegas vs x
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(x_array, OmegaR, label=r'$\Omega_R$')
    ax.plot(x_array, OmegaCDM, label=r'$\Omega_{CDM}$')
    ax.plot(x_array, OmegaB, label=r'$\Omega_B$')
    ax.plot(x_array, OmegaLambda, label=r'$\Omega_\Lambda$')
    ax.plot(x_array, Omega_matter, label=r'$\Omega_{\rm{Matter}}$', ls='--', alpha=0.7)
    color_each_regime(ax,x_array)
    ax.margins(x=0)
    ax.legend(loc=6)
    ax.set_xlabel(r'$x = ln(a)$')
    ax.set_ylabel(r'$\Omega_i = \frac{\rho_i}{\rho_c}$')
    ax.set_xticks(x_array_ticks)

    print('Saving ./../figs/omegas_of_x.pdf')
    fig.savefig('./../figs/omegas_of_x.pdf', bbox_inches='tight')

    # Plotting Hubble parameter vs x
    fig, axes = plt.subplots(2, 2,figsize=(12,9))
    axes[0, 0].plot(x_array, H_units, color='C2')
    axes[0, 0].set_ylabel(r'$H(x) \left[\rm{km}\,\rm{s^{-1}}\,\rm{Mpc^{-1}}\right]$')

    # Plotting Hubble parameter vs z
    axes[0, 1].plot(z_non_neg, H_units[z_array > 0], color='C2')
    axes[0, 1].set_xlabel(r'$z = a^{-1} - 1$')
    axes[0, 1].set_ylabel(r'$H(z) \left[\rm{km}\,\rm{s^{-1}}\,\rm{Mpc^{-1}}\right]$')
    axes[0, 1].set_xlim(z_non_neg.max(), z_non_neg.min()+1e-8)
    axes[0, 1].set_xscale('log')
    axes[0, 1].set_xticks(z_array_ticks)
    color_each_regime(axes[0,1],x_array,z_array)

    # Plotting conformal time
    axes[1, 0].plot(x_array, eta_units, color='C2')
    axes[1, 0].plot(x_array, aeta_R, color='C0', ls='-.')
    axes[1, 0].plot(x_array, aeta_M, color='C1', ls='-.')
    axes[1, 0].plot(x_array, aeta_Lambda, color='C3', ls='-.')
    axes[1, 0].set_ylabel(r'$\eta(x) [\rm{Gpc}]$')

    # PLotting H-prime
    axes[1, 1].plot(x_array, Hp_units, color='C2')
    axes[1, 1].set_ylabel(r'$\mathcal{H}(x) \left[\rm{km}\,\rm{s^{-1}}\,\rm{Mpc^{-1}}\right]$')

    # Filling background in each axis
    for i, ax in enumerate(axes.flatten()):
        ax.set_yscale('log')
        ax.margins(x=0)
        if i != 1:
            color_each_regime(ax,x_array)
            ax.set_xlabel(r'$x = ln(a)$')
            ax.set_xticks(x_array_ticks)
    print('Saving ./../figs/Hubble_eta_of_x.pdf')
    fig.savefig('./../figs/Hubble_eta_of_x.pdf', bbox_inches='tight')
    

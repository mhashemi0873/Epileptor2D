def calculate_summary_statistics(x, dt, ts):
    """Calculate summary statistics

    Parameters
    ----------
    x : output of the simulator

    Returns
    -------
    np.array, summary statistics
    """


    t_on=ts[0]
    t_off=ts[-1]

    n_mom = 4
    n_summary=7
    n_summary = np.minimum(n_summary, n_mom + 3)

    # initialise array of spike counts
    v=np.zeros(nt)
    v= np.array(x)


    # put everything to 0 that is below 0 or has negative slope
    ind = np.where(v < 0)
    v[ind] = 0
    ind = np.where(np.diff(v) < 0)
    v[ind] = 0

    # remaining negative slopes are at spike peaks
    ind = np.where(np.diff(v) < 0)

    spike_times = np.array(ts)[ind]
    spike_times_stim = spike_times[(spike_times > t_on) & (spike_times < t_off)]


    if spike_times_stim.shape[0] > 0:
        spike_times_stim = spike_times_stim[np.append(1, np.diff(spike_times_stim)) > 0.5]

        
    std_pw = np.power(np.std(x[(ts > t_on) & (ts < t_off)]), np.linspace(3, n_mom, n_mom - 2) )

    std_pw = np.concatenate((np.ones(1), std_pw))


    #momments
    spstats.moment(x[(ts > t_on) & (ts < t_off)], np.linspace(2, n_mom, n_mom - 1))/std_pw

    moments = spstats.moment(v, np.linspace(2, n_mom, n_mom - 1))/ std_pw

    sum_stats_vec = np.concatenate((np.array([spike_times_stim.shape[0]]),
                                    np.array([np.mean(x)]),
                                    moments,))

    sum_stats_vec = sum_stats_vec[0:n_summary]
    


    return sum_stats_vec

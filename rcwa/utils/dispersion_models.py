def drudeLorentzModel(wavelength, plasma_frequency, damping_coefficients,
        oscillator_amplitudes, offset_coefficients):
    """
    :param wavelength: The frequency of interest (in microns)
    :param plasma_frequency: The plasma frequency of the metal of interest
    :param damping_coefficients: A k+1 array of damping coefficients ( in electronvolts), k is number of oscillators.
    :param oscillator_amplitudes: A k+1 array of oscillator amplitudes (unitless). Zeroth element taken to be amplitude of intraband transitions.
    :param offset_coefficients: A k array of offset frequencies (in electronvolts)
    :returns: A two-tuple containing permittivity and permeability (er, ur)
    """

    frequency = 1.2398419843320028 / wavelength # converts from microns to eV
    er_free = 1 - oscillator_amplitudes[0] * plasma_frequency * plasma_frequency / \
            (frequency * (frequency - 1j*damping_coefficients[0]))
    er_bound = 0
    for i in range(len(damping_coefficients - 1)):
        er_bound += oscillator_amplitudes[i] * plasma_frequency * plasma_frequency / \
                (offset_coefficients*offset_coefficients - \
                frequency*frequency + 1j*frequency*damping_coefficients[i])
    er = er_free + er_bound
    ur = 1
    return (er, ur)

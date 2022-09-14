import astropy.units as u

def freq_rest():

    ### NIST https://www.nist.gov/pml/observed-interstellar-molecular-microwave-transitions
    ### 

    freq_rest = {'12CO_1-0': 115.271202 * u.GHz,  ## (v=0) NIST / ALMA OT
                 '12CO_2-1': 230.538000 * u.GHz,  ## (v=0) NIST / ALMA OT
                 '12CO_3-2': 345.795990 * u.GHz,  ## (v=0) NIST / ALMA OT
                 '12CO_4-3': 461.040768 * u.GHz,  ## (v=0) NIST / ALMA OT
                 '12CO_5-4': 576.267931 * u.GHz,  ## (v=0) NIST / ALMA OT
                 '12CO_6-5': 691.473076 * u.GHz,  ## (v=0) NIST / ALMA OT
                 '12CO_7-6': 806.651801 * u.GHz,  ## (v=0) NIST / ALMA OT
                 '12CO_8-7': 921.799704 * u.GHz,  ## (v=0) NIST / ALMA OT
                 '12CO_9-8': 1036.912385 * u.GHz, ## (v=0) NIST
                 '12CO_10-9': 1151.985 * u.GHz,   ## (v=0) Alcolea et al. (2013) https://ui.adsabs.harvard.edu/abs/2013A%26A...559A..93A/abstract
                 '12CO_11-10': 1267.014482 * u.GHz,  ## (v=0) NIST
                 '12CO_12-11': 1381.995102 * u.GHz,  ## (v=0) NIST
                 '12CO_14-13': 1611.793508 * u.GHz,  ## (v=0) NIST
                 '12CO_11-10': 1267.014482 * u.GHz,  ## (v=0) NIST 
                 '12CO_16-15': 1841.345512 * u.GHz,  ## (v=0) NIST / Alcolea et al. (2013) https://ui.adsabs.harvard.edu/abs/2013A%26A...559A..93A/abstract
                 '12CO_17-16': 1956.018137 * u.GHz,  ## (v=0) NIST / Gallerani et al. (2014) https://ui.adsabs.harvard.edu/abs/2014MNRAS.445.2848G/abstract
                 }

        
    return freq_rest

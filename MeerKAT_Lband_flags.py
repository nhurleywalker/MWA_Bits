def apply_flags(data):
    data[:,37:86] = np.nan
    data[:,165:170] = np.nan
    data[:,229:269] = np.nan
    data[:,287:506] = np.nan
    data[:,585:590] = np.nan
    data[:,747:856] = np.nan
    data[:,868:884] = np.nan
    data[:,921:] = np.nan
    return(data)


def print_flags():
    starts = [37, 165, 229, 287, 585, 747, 868, 921]
    ends = [86, 170, 269, 506, 590, 856, 884, 932]
    for s, e in zip(starts, ends):
        print(np.arange(s, e))


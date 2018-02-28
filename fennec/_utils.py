
def _max_proc(min_cpu=1):
    '''
    Return the maximum number of available CPU defined as:
        number of cpu - roundup(average load over the last minute)

    Will return at least `min_cpu` (1 by default)
    '''
    import os
    return max(os.cpu_count() - int(os.getloadavg()[0] + 1), min_cpu)


def _print_progressbar(step, maxi, msg="", char="=", width=50):
    '''
    Print a progress bar then place the cursor at the begging of the line.
    Display can be really messy if `maxi` is set incorrectly.

    import timeg    n=32
    for i in range(n):
        time.sleep(0.1)
        _print_progressbar(i+1, n, msg="Test", char='=', width=50)
    print()

    '''
    # rows, columns = os.popen('stty size', 'r').read().split()
    p = int(100 * step / maxi)
    print("[%s>%s] %d%% (%d/%d) %-20s" %
        (char * int(p * width/100),
        (" " * (width-int(p*width/100))),
        p, step, maxi, msg),
        end="\r", flush=True)

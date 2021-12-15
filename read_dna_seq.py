def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return numpy.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
handle = open("U00096.2.fas.txt", "r")
a = handle.readlines()[1:]
a = ''.join([x.strip() for x in a])
conv = {'a': 0, 'c': 1, 'g': 2, 't': 3}
b = numpy.array([conv[x] for x in a], dtype=numpy.uint8)
rolling_window(b,100)

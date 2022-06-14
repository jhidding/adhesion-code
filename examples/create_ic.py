# ~\~ language=Python filename=examples/create_ic.py
# ~\~ begin <<tutorial.md|create-ic>>[init]
import numpy as np
import h5py as h5
# ~\~ end
# ~\~ begin <<tutorial.md|create-ic>>[1]
L = 50.0
N = 128
shape = (N, N, N)

white_noise = np.random.normal(0.0, 1.0, shape)
# ~\~ end
# ~\~ begin <<tutorial.md|create-ic>>[2]
def _wave_number(s):
    N = s[0]
    i = np.indices(s)
    return np.where(i > N/2, i - N, i)

K = _wave_number(shape) * 2*np.pi / L
k = np.sqrt((K**2).sum(axis=0))
# ~\~ end
# ~\~ begin <<tutorial.md|create-ic>>[3]
np_default_err = np.seterr(divide = 'ignore', invalid = 'ignore')
density_fourier = np.fft.fftn(white_noise) * k**-1.5
density_fourier.flat[0] = 0
density = np.fft.ifftn(density_fourier).real
density /= density.std()
# ~\~ end
# ~\~ begin <<tutorial.md|create-ic>>[4]
potential_fourier = np.fft.fftn(density) * k**-2
potential_fourier.flat[0] = 0
potential = np.fft.ifftn(potential_fourier).real
np.seterr(**np_default_err)
# ~\~ end
# ~\~ begin <<tutorial.md|create-ic>>[5]
with h5.File("my-ic.h5", "w") as f:
    f.attrs["N"] = N
    f.attrs["L"] = L
    f.create_dataset("potential", data=potential)
# ~\~ end

# Tutorial: use your own potential
You may want to try the adhesion model on your own initial conditions. In that case you need to prepare an HDF5 file that contains the attributes for physical size `L`, the logical size `N` and a data set of shape `(N,N,N)` named `potential`. The easiest way to create such a file is in Python.

## Dependencies
Create an environment (using Conda or Poetry) that has both NumPy and H5Py installed.

``` {.python file=examples/create_ic.py #create-ic}
import numpy as np
import h5py as h5
```

## Generate random Gaussian noise
We will generate some scale-free density perturbations with a spectrum of $P(k) \sim k^{-3}$. Start by making some white noise.

``` {.python #create-ic}
L = 50.0
N = 128
shape = (N, N, N)

white_noise = np.random.normal(0.0, 1.0, shape)
```

To apply a power spectrum, we need to compute the values of the Fourier frequencies on our grid.

``` {.python #create-ic}
def _wave_number(s):
    N = s[0]
    i = np.indices(s)
    return np.where(i > N/2, i - N, i)

K = _wave_number(shape) * 2*np.pi / L
k = np.sqrt((K**2).sum(axis=0))
```

Now we can compute the initial density field. Take care not to divide by zero. We normalize this field by setting the standard deviation to 1. This corresponds roughly to having our first structures collapse at time D=1.

``` {.python #create-ic}
np_default_err = np.seterr(divide = 'ignore', invalid = 'ignore')
density_fourier = np.fft.fftn(white_noise) * k**-1.5
density_fourier.flat[0] = 0
density = np.fft.ifftn(density_fourier).real
density /= density.std()
```

Now we should compute the potential.

``` {.python #create-ic}
potential_fourier = np.fft.fftn(density) * k**-2
potential_fourier.flat[0] = 0
potential = np.fft.ifftn(potential_fourier).real
np.seterr(**np_default_err)
```

## Write output and run
Write this to HDF5

``` {.python #create-ic}
with h5.File("my-ic.h5", "w") as f:
    f.attrs["N"] = N
    f.attrs["L"] = L
    f.create_dataset("potential", data=potential)
```

To run the adhesion model on our custom ICs, we need to write a configuration file:

``` {.yaml file=examples/ic-from-file.yaml}
initial-conditions:
  file: my-ic.h5

run:
  time:     [0.5, 1.0, 2.0]

output:
  hdf5:            ./my-ic.h5
  walls:           output/my-ic-{time:02.1f}-walls.obj
  filaments:       output/my-ic-{time:02.1f}-filaments.obj
  threshold:       1.0
```

And now we can run the adhesion code (in this case from the examples directory). Make sure that the `output` directory exists, otherwise nothing will get written (should fix):

```bash
python ./create-ic.py
mkdir output
../build/adhesion -c ic-from-file.yaml
```


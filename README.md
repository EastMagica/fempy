# fempy
A simple Python-based Finite Element Method Library.

## Schedule

* 2019.07.12 getting start
* 2019.07.13 Finished a simple 1-D Finity Element Method 

## Getting start

### prerequisites

* python 3.6
* numpy 1.16.3
* scipy 1.0.1
* matplotlib 3.0.2

### installing

Installation is simple! clone this repo and run

```shell
python setup.py install
```

at the root folder.

## Example

### 1-D simple equation

$$\begin{cases}
-u''=\pi^2\sin(\pi x) \\
u(0)=0, u(1)=0
\end{cases}$$

we konw that the true solution of u is $u(x)=\sin(\pi x)$, the using: 

```python
eq = {'bnd': np.array([0, 0]),
      'f': lambda x: np.pi**2 * np.sin(np.pi * x),
      'domain': np.array([0, 1]), 'split': ['step', 8]}
fem = FEM1D(eq)
fem.run()
```

we could know the numerical solution obtained by finite element method.

## Authors

* **EastMagica** - *Initial work* - [EastMagica](https://github.com/EastMagica)

## License
This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details

## Acknowledgments
* Yi Nianyu, Xiangtan University
* Wang dongling, Northwest University
* etc

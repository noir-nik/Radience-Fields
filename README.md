# Radiance Fields with Spherical Harmonics

This project implements a method of reconstructing 3D scenes from a set of images described in [Plenoxels: Radiance Fields without Neural Networks](https://arxiv.org/abs/2112.05131) and [ReLU Fields: The Little Non-linearity That Could](https://arxiv.org/abs/2205.10824)

## Features
 - Rendering images from a voxel grid where each voxel stores an opacity and spherical harmonics coefficients for each color channel.
 - Training a radiance field model from a dataset of images

## Rendering

![](/images/rendering_model.jpg)

### Building render target

```bash
git clone https://github.com/noir-nik/Radience-Fields.git
cd render
make -j8
```

### Usage
```bash
rf-render -m <path to binary grid> -s <grid size>
```
Grid is stored as a dense array of `struct Cell` of size `grid_size * grid_size * grid_size`
```c
struct Cell
{
	float density;
	float sh_r[SH_WIDTH];
	float sh_g[SH_WIDTH];
	float sh_b[SH_WIDTH];
};
```

`SH_WIDTH` is defined in `settings.h`

The example model can be used with the following command:

```bash
rf-render -m data/model.dat -s 128
```

The images will be saved in the `output` folder

## Training
![](/images/traning_process.jpg)
### Building training target

The training target requires clang compiler and the [Enzyme AD](https://github.com/EnzymeAD/Enzyme) library to be installed.

```bash
cd training
make
```
### Running:
```bash
make run
```
Choose one of the available datasets: type `lego` or `chair`

Training is done with a coarse-to-fine strategy, starting with a grid size of 64 and then increasing to 128.

The number of references from the dataset (from 1 to 100) can be changed in the `settings.h` file, as well as other training parameters.

It is recommended to have more than 8 GB of free RAM.

When the model converges, the grid will be saved to the 'output_grid_data' folder.
Intermediate images and the final render after training will appear in the `output` folder.
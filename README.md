# CS171 Final Project —— Ray Tracing NURBS surface

Greetings from Group 10 (Xiyue Peng `<pengxy1>` and Luojia Hu `<hulj>`) !

## Usage

### Libs

We used `Eigen`, `nlohmann::json`, `stb` and `tinyobjloader` as our libs. Those libs have to be placed inside `libs/` folder in order to compile the project.

### Executable

We provided an executable in `build` folder. It was compiled it on an Intel-chip MacBook Pro with gcc-12 on macOS Ventura, so it should be able to run on x86-64 macOS devices. Make sure you are in the correct directory when execute it, otherwise it may fail to read the config files (located in `configs/`) and object parameters (located in `assets/`).

We have provided a NURBS sphere in `assets/sphere.nurbs`. We also set up a scene with the NURBS object inside it. Run the following to get the rendered image.

```shell
cd build
./CS171-final-project ../configs/final_project.json
```

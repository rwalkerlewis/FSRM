# FSRM Development Environment

This directory contains configuration for VS Code Dev Containers, providing a consistent, cross-platform development environment for FSRM.

## Quick Start

### Option 1: VS Code Dev Containers (Recommended)

1. **Install Prerequisites:**
   - [Docker Desktop](https://www.docker.com/products/docker-desktop/)
   - [VS Code](https://code.visualstudio.com/)
   - [Dev Containers Extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

2. **Open in Container:**
   - Open the FSRM folder in VS Code
   - Click the green button in the bottom-left corner (or press `F1`)
   - Select "Dev Containers: Reopen in Container"
   - Wait for the container to build (first time takes ~10-15 minutes)

3. **Build and Test:**
   ```bash
   cd /workspace/build
   make -j$(nproc)
   ctest --output-on-failure
   ```

### Option 2: Docker Compose

```bash
# Build and start the development container
docker compose -f docker-compose.dev.yml up -d

# Enter the container
docker compose -f docker-compose.dev.yml exec fsrm-dev bash

# Build FSRM inside the container
cd /workspace/build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_TESTING=ON
make -j$(nproc)
ctest
```

### Option 3: Manual Docker

```bash
# Build the development image
docker build -f Dockerfile.dev -t fsrm-dev:latest .

# Run the container
docker run -it --rm \
  -v $(pwd):/workspace \
  --cap-add=SYS_PTRACE \
  --security-opt seccomp=unconfined \
  --ipc=host \
  fsrm-dev:latest bash

# Build inside container
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_TESTING=ON
make -j$(nproc)
```

## What's Included

### Development Tools
- **Build System:** CMake, Ninja, Make
- **Compilers:** GCC, GFortran
- **Debugger:** GDB, Valgrind
- **Version Control:** Git, Git LFS

### Libraries
- **PETSc:** Parallel solver framework
- **MPI:** OpenMPI for parallel computing
- **HDF5:** Data I/O with MPI support
- **PROJ:** Coordinate transformations
- **Gmsh:** Mesh generation

### Python Environment
- NumPy, SciPy, Matplotlib
- h5py, pandas
- pytest for testing
- JupyterLab for notebooks

### VS Code Extensions (auto-installed)
- C/C++ IntelliSense
- CMake Tools
- GitLens
- Debugger support
- GitHub Copilot

## Useful Commands

Inside the container, these aliases are available:

```bash
# Quick build (from any directory)
build

# Run tests
test

# Reconfigure CMake
configure
```

## Debugging

The container is configured for debugging with GDB:

1. Set breakpoints in VS Code
2. Press `F5` or use "Run and Debug"
3. Select "C++ (GDB)" configuration

For MPI debugging:
```bash
mpirun -np 4 xterm -e gdb ./my_program
```

## Customization

### Adding Dependencies

Edit `Dockerfile.dev` to add system packages:
```dockerfile
RUN apt-get update && apt-get install -y \
    your-package-here
```

### VS Code Settings

Modify `.devcontainer/devcontainer.json` to add extensions or settings.

### Environment Variables

Add to `docker-compose.dev.yml` or `devcontainer.json`:
```yaml
environment:
  - MY_VAR=value
```

## Troubleshooting

### Container won't start
```bash
# Check Docker is running
docker info

# Remove old containers and rebuild
docker compose -f docker-compose.dev.yml down
docker compose -f docker-compose.dev.yml build --no-cache
```

### Build errors
```bash
# Clean build directory
rm -rf /workspace/build/*
cmake .. -DCMAKE_BUILD_TYPE=Debug
```

### MPI errors
```bash
# Ensure IPC is shared
docker run --ipc=host ...
```

### Permission issues
The container runs as user `developer` (UID 1000). If you have permission issues:
```bash
# Inside container
sudo chown -R developer:developer /workspace
```

# FSRM Build Script for Docker on Windows
# Prerequisites: Docker Desktop must be installed and running
# Download from: https://www.docker.com/products/docker-desktop/

Write-Host "=================================" -ForegroundColor Cyan
Write-Host "FSRM Build Script for Docker" -ForegroundColor Cyan
Write-Host "=================================" -ForegroundColor Cyan
Write-Host ""

# Check if Docker is available
try {
    $dockerVersion = docker --version 2>&1
    Write-Host "Found: $dockerVersion" -ForegroundColor Green
} catch {
    Write-Host "ERROR: Docker is not installed or not in PATH" -ForegroundColor Red
    Write-Host "Please install Docker Desktop from: https://www.docker.com/products/docker-desktop/" -ForegroundColor Yellow
    exit 1
}

# Navigate to project directory
Set-Location "c:\Users\User\Projects\FSRM"

Write-Host ""
Write-Host "Step 1: Building Docker image (this may take 10-20 minutes)..." -ForegroundColor Cyan
docker build -t fsrm:latest .

if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Docker build failed" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "Step 2: Creating output directories..." -ForegroundColor Cyan
New-Item -ItemType Directory -Force -Path ".\output" | Out-Null
New-Item -ItemType Directory -Force -Path ".\data" | Out-Null
New-Item -ItemType Directory -Force -Path ".\simulations" | Out-Null

Write-Host ""
Write-Host "Step 3: Starting container..." -ForegroundColor Cyan
docker-compose up -d

if ($LASTEXITCODE -ne 0) {
    Write-Host "WARNING: docker-compose failed, trying direct docker run..." -ForegroundColor Yellow
    docker run -d --name fsrm `
        -v "${PWD}/config:/config" `
        -v "${PWD}/output:/app/output" `
        -v "${PWD}/data:/data" `
        fsrm:latest tail -f /dev/null
}

Write-Host ""
Write-Host "Step 4: Running tests inside container..." -ForegroundColor Cyan
docker exec fsrm bash -c "cd /app/build && ctest --output-on-failure"

Write-Host ""
Write-Host "=================================" -ForegroundColor Green
Write-Host "Build Complete!" -ForegroundColor Green
Write-Host "=================================" -ForegroundColor Green
Write-Host ""
Write-Host "Container 'fsrm' is now running." -ForegroundColor White
Write-Host ""
Write-Host "To run examples:" -ForegroundColor Yellow
Write-Host "  docker exec -it fsrm bash" -ForegroundColor White
Write-Host "  cd /app/build/examples" -ForegroundColor White
Write-Host "  ./simulator -c /app/config/default.config" -ForegroundColor White
Write-Host ""
Write-Host "To run with MPI:" -ForegroundColor Yellow
Write-Host "  docker exec fsrm mpirun -np 4 /app/build/examples/simulator -c /app/config/hydraulic_fracturing.config" -ForegroundColor White
Write-Host ""
Write-Host "To stop the container:" -ForegroundColor Yellow
Write-Host "  docker stop fsrm" -ForegroundColor White
Write-Host ""
Write-Host "Output files will be in: .\output\" -ForegroundColor Yellow
Write-Host ""

@echo off
echo Starting CompuCell3D mitosis_O2 simulation with VTK output...
echo.

REM Clean up old output directory INSIDE Cell_Spheroid working directory
if exist "output" (
    echo Removing old output folder...
    rmdir /s /q "output"
    echo Cleaned up old results
)

REM Create the output directory structure INSIDE Cell_Spheroid
mkdir "output"

REM Run the simulation with output directory INSIDE Cell_Spheroid working directory
echo Running simulation...
C:\CompuCell3D\runScript.bat -i mitosis_O2.cc3d -o output -f 50 --current-dir="%CD%" > simulation.log 2>&1

echo.
if %errorlevel% equ 0 (
    echo Simulation completed successfully!
    
    REM Count VTK files created IN the Cell_Spheroid working directory
    if exist "output\Step_*.vtk" (
        for /f %%i in ('dir /b output\Step_*.vtk 2^>nul ^| find /c /v ""') do set vtk_count=%%i
        echo Created %vtk_count% VTK files in .\output\
    ) else (
        echo Warning: No VTK files found in .\output\
    )
) else (
    echo Simulation failed with error code %errorlevel%
)

echo.
echo Log file: simulation.log
echo VTK files are saved LOCALLY in: .\output\
echo Output is INSIDE your Cell_Spheroid working directory!
echo You can now open the VTK files in ParaView.
echo.
echo ParaView Instructions:
echo 1. Open ParaView
echo 2. File → Open → Navigate to .\output\ folder (INSIDE Cell_Spheroid)
echo 3. Select all Step_*.vtk files (they will group automatically)
echo 4. Click Apply
echo 5. Change coloring to 'CellType' to see the 3 phenotypes
echo 6. Use the play button to animate through time
echo 7. For 3D visualization: Change representation to 'Volume' or use 'Contour' filter
echo.
pause
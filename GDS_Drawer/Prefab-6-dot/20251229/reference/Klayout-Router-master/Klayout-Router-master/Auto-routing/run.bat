@echo off
REM ==========================================
REM  Auto-Klayout Routing Configuration Script
REM ==========================================

REM --- File and Cell Settings ---
set FILE_PATH="D:/_Topological quantum computation/PbTe/Quantum Point Contact (QPC)/SC487 - QPC above mask/SC487 - 副本.gds"
set CELL_NAME="Top"

REM --- Layer Settings ---
set BBOX_LAYER="7/0"
set PIN_A_LAYER="102/1"
set PIN_B_LAYER="111/1"
set OUTPUT_LAYER="121/0"

REM --- Map Settings ---
set RESOLUTION=5

REM --- Routing Parameters ---
set OBS_SAFE_DISTANCE=16
set OBS_HARDNESS=20
set OBS_DAMPING_STEP=4

set PIN_A_SAFE_DIST=20
set PIN_B_SAFE_DIST=200
set PIN_HARDNESS=20
set PIN_A_DAMPING=4
set PIN_B_DAMPING=1

set WIDTH=4
set PATH_SAFE_DIST=20
set PATH_HARDNESS=10
set PATH_DAMPING=5
set PATH_DENSITY_HARDNESS=1

set ROUND=1

REM --- Animation Settings ---
REM Set ANIMATE_FLAG to "--animate" to enable, or leave blank "" to disable
set ANIMATE_FLAG=--animate
set SAVE_ANIM="routing.gif"
set FPS=20

REM --- Execution ---
echo Running Auto-Klayout Routing...
python routing.py ^
    --file_path %FILE_PATH% ^
    --cell_name %CELL_NAME% ^
    --bbox_layer %BBOX_LAYER% ^
    --pin_a_layer %PIN_A_LAYER% ^
    --pin_b_layer %PIN_B_LAYER% ^
    --output_layer %OUTPUT_LAYER% ^
    --resolution %RESOLUTION% ^
    --obs_safe_distance %OBS_SAFE_DISTANCE% ^
    --obs_hardness %OBS_HARDNESS% ^
    --obs_damping_step %OBS_DAMPING_STEP% ^
    --pin_a_safe_dist %PIN_A_SAFE_DIST% ^
    --pin_b_safe_dist %PIN_B_SAFE_DIST% ^
    --pin_hardness %PIN_HARDNESS% ^
    --pin_a_damping %PIN_A_DAMPING% ^
    --pin_b_damping %PIN_B_DAMPING% ^
    --width %WIDTH% ^
    --path_safe_dist %PATH_SAFE_DIST% ^
    --path_hardness %PATH_HARDNESS% ^
    --path_damping %PATH_DAMPING% ^
    --path_density_hardness %PATH_DENSITY_HARDNESS% ^
    --round %ROUND% ^
    %ANIMATE_FLAG% ^
    --save_anim %SAVE_ANIM% ^
    --fps %FPS%

if %ERRORLEVEL% NEQ 0 (
    echo Error occurred during execution.
    pause
    exit /b %ERRORLEVEL%
)

echo Done.
pause

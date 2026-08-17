@echo off
chcp 65001 >nul
cd /d "%~dp0"

echo ==============================================
echo   多免疫检查点小分子AI筛选系统 - 启动器
echo   使用 conda 环境: immuno_drug_screen
echo ==============================================
echo.

set CONDA_ENV_PATH=%USERPROFILE%\miniconda3\envs\immuno_drug_screen

REM 检查 conda 环境是否存在
if not exist "%CONDA_ENV_PATH%\python.exe" (
    echo [错误] 未找到 conda 环境 immuno_drug_screen
    echo   请先运行: conda create -n immuno_drug_screen python=3.13 pip -y
    echo   然后安装依赖: pip install flask numpy pandas scipy pyyaml scikit-learn rdkit openbabel torch torch-geometric e3nn ema-pytorch fair-esm torchmetrics prody
    pause
    exit /b 1
)

echo [1/2] 使用 conda 环境 Python: %CONDA_ENV_PATH%\python.exe
echo [2/2] 启动 Web 服务器...
echo.
echo   本机访问:   http://127.0.0.1:5050
echo   局域网访问: http://您的IP:5050
echo   按 Ctrl+C 停止服务器
echo.
echo ==============================================
echo.

"%CONDA_ENV_PATH%\python.exe" app.py

pause
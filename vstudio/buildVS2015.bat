call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\amd64\vcvars64.bat"

"C:\Program Files (x86)\MSBuild\14.0\Bin\MSBuild.exe" Splat.sln /p:configuration=Debug /p:platform=x64 /t:clean
"C:\Program Files (x86)\MSBuild\14.0\Bin\MSBuild.exe" Splat.sln /p:configuration=Release /p:platform=x64 /t:clean
"C:\Program Files (x86)\MSBuild\14.0\Bin\MSBuild.exe" Splat.sln /p:configuration=Release /p:platform=x64


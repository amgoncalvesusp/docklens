; DockLens Windows installer.
; The CI build passes absolute SourceExe and OutputDir values so the installer
; wraps the already smoke-tested PyInstaller executable.

#ifndef AppVersion
  #define AppVersion "1.2.0"
#endif
#ifndef SourceExe
  #define SourceExe "..\dist\DockLens.exe"
#endif
#ifndef OutputDir
  #define OutputDir "..\dist"
#endif

[Setup]
AppId={{B7F1C6D0-4F92-4E31-9B0D-4E95C2F16A24}
AppName=DockLens
AppVersion={#AppVersion}
AppPublisher=Adriano Marques Gonçalves — UNIARA
AppPublisherURL=https://github.com/amgoncalvesusp/docklens
AppSupportURL=https://github.com/amgoncalvesusp/docklens/issues
DefaultDirName={autopf}\DockLens
DefaultGroupName=DockLens
DisableProgramGroupPage=yes
OutputDir={#OutputDir}
OutputBaseFilename=DockLens-setup-windows-x86_64
Compression=lzma2
SolidCompression=yes
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
PrivilegesRequired=admin
UninstallDisplayName=DockLens
WizardStyle=modern

[Files]
Source: "{#SourceExe}"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{autoprograms}\DockLens"; Filename: "{app}\DockLens.exe"
Name: "{autodesktop}\DockLens"; Filename: "{app}\DockLens.exe"

[Run]
Filename: "{app}\DockLens.exe"; Description: "Launch DockLens"; Flags: postinstall nowait skipifsilent

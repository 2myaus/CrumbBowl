{ pkgs ? import <nixpkgs> {} }:
  pkgs.mkShell {
    # nativeBuildInputs is usually what you want -- tools you need to run
    nativeBuildInputs = with pkgs.buildPackages; [
    libgcc
    # libclang
    clang-tools
    gdb
    valgrind
    kdePackages.kcachegrind
    ];
}

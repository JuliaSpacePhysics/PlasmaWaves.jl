# Install ALPS
install_ALPS:
    #!/usr/bin/env sh
    [ -d ALPS ] || git clone https://github.com/danielver02/ALPS.git
    cd ALPS
    ./configure
    make

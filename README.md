Process:
1. setting the port in docker-compose.yaml
2. docker-compose up --build
3. cd /root/projects/libRadtran-2.0.4
>3. download the other version
    1. download libRadtran from offical website: http://www.libradtran.org/doku.php?id=download
    2. gzip libRadtran-{version}.tar.gz
    3. tar -vxf libRadtran-{version}.tar
    4. cd {path of libRadtran}
7. ./configure
8. make
9. make check
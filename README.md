Process:
1. setting the port in docker-compose.yaml
2. docker-compose up --build
3. download libRadtran from offical website: http://www.libradtran.org/doku.php?id=download
4. gzip libRadtran-{version}.tar.gz
5. tar -vxf libRadtran-{version}.tar
6. cd {path of libRadtran}
7. ./configure
8. make
9. make check
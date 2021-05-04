Process:
1. setting the port in docker-compose.yaml
2. docker-compose up --build
3. gzip libRadtran-{version}.tar.gz
4. tar -vxf libRadtran-{version}.tar
5. cd {path of libRadtran}
6. ./configure
7. make
8. make check
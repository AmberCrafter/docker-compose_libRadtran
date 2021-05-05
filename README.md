Process:
1. setting the port in docker-compose.yaml
2. setting the username and password in ubuntu20/dockerfile #4 and #5
3. docker-compose up --build
4. cd /root/projects/libRadtran-2.0.4
  >4. download the other version
    1. download libRadtran from offical website: http://www.libradtran.org/doku.php?id=download
    2. gzip libRadtran-{version}.tar.gz
    3. tar -vxf libRadtran-{version}.tar
    4. cd {path of libRadtran}
5. ./configure
6. make
7. make check
  >5. configure and make by shell script (with check)
    1. cd ..
    2. bash build.sh
FROM alpine

# libtbb-dev is in edge/testing, and can't be properly installed by vcpkg
RUN echo "https://dl-cdn.alpinelinux.org/alpine/edge/testing" >> /etc/apk/repositories
RUN apk update
RUN apk add make g++ cmake sqlite-dev libtbb-dev

RUN mkdir /Isconna
WORKDIR /Isconna
COPY CMakeLists.txt .
COPY dep dep
COPY include include
COPY example example

RUN cmake -DCMAKE_BUILD_TYPE=Release -S . -B build/release
RUN cmake --build build/release --target Demo

ENTRYPOINT ["build/release/Demo"]

FROM alpine/git

# libtbb-dev is in edge/testing, and can't be properly installed by vcpkg
RUN echo "https://dl-cdn.alpinelinux.org/alpine/edge/testing" >> /etc/apk/repositories
RUN apk update
RUN apk add make g++ cmake sqlite-dev libtbb-dev

RUN mkdir /Isconna
WORKDIR /Isconna
RUN git clone https://github.com/liurui39660/AUROC dep/AUROC --depth=1
RUN git clone https://github.com/mandreyel/mio.git dep/mio --depth=1
COPY example example
COPY src src
COPY CMakeLists.txt .

RUN cmake -DCMAKE_BUILD_TYPE=Release -S . -B build/release
RUN cmake --build build/release --target Demo

ENTRYPOINT ["build/release/Demo"]

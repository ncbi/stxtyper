#!/bin/bash
# found and lightly modified these functions
get_latest_release() {
    curl --silent "https://api.github.com/repos/$1/releases/latest" |
    grep '"tag_name":' |
    cut -d '"' -f 4
}

get_tarball_url() {
    curl --silent "https://api.github.com/repos/$1/releases/latest" |
        fgrep '"browser_download_url":' |
        cut -d '"' -f 4
}


release=$(get_latest_release ncbi/stxtyper)
URL=$(get_tarball_url ncbi/stxtyper)

>&2 echo "Downloading StxTyper version $release"
>&2 echo "Binaries from $URL"

# download and unpack AMRFinder binaries
    curl --silent -L -O $URL
    tarball_name=$(echo $URL | perl -pe 's#^.*/(.*)#\1#')
    tar xfz $tarball_name
    rm $tarball_name
    directory=$(echo $tarball_name | sed 's/.tar.gz//')
    echo $directory

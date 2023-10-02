
set -e
set -o pipefail

# Use the OpenPBTA bucket as the default.
URL=${PBTA_URL:-https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/pbta-ancestry}
RELEASE=${PBTA_RELEASE:-v1}

# Remove old symlinks in data
find data -type l -delete

# The md5sum file provides our single point of truth for which files are in a release.
curl --create-dirs $URL/$RELEASE/md5sum.txt -o data/$RELEASE/md5sum.txt -z data/$RELEASE/md5sum.txt

# Consider the filenames in the md5sum file and the release notes
FILES=(`tr -s " " < data/$RELEASE/md5sum.txt | cut -d$'\t' -f 2` release-notes.md)

# Download the items in FILES if not already present
for file in "${FILES[@]}"
do
  if [ ! -e "data/$RELEASE/$file" ]
  then
    echo "Downloading $file"
    curl $URL/$RELEASE/$file -o data/$RELEASE/$file
  fi
done


# Check the md5s for everything we downloaded except CHANGELOG.md
cd data/$RELEASE
echo "Checking MD5 hashes..."
md5sum -c md5sum.txt
cd ../../

# Make symlinks in data/ to the files in the just downloaded release folder.
for file in "${FILES[@]}"
do
  ln -sfn $RELEASE/$file data/$file
done

# make data directory unwritable in CI
if [ "$RELEASE" == "testing" ]; then
  chmod u-w data
fi

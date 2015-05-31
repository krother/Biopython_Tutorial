
# example phrases useful for building pipelines:

# do something for all files in a directory
import os

for filename in os.listdir("."):
    if filename.endswith(".py"):
        print filename


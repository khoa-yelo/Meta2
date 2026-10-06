# checkgm scoring image: the Python package only (no assembler or profiler). Use it to normalize and score MGnify v5 analyses,
# query directories or MetaPhlAn 3 / HUMAnN 3 tables against a bundle mounted into the container.
#   podman build -t checkgm:0.8.2 -f container/Dockerfile .          # run from the checkgm/ tool directory
#   (on a rootless host without subuid ranges add: --cgroup-manager=cgroupfs --events-backend=file)
#   podman run --rm -v /data/checkgm_bundles:/bundles:ro -v $PWD:/work -e CHECKGM_BUNDLES=/bundles checkgm:0.8.2 \
#       run-metaphlan --bundle gut-reads-adult-global-v0.2-lenient --input /work/S1.txt --reads-tsv /work/reads.tsv --out /work/S1_report
FROM docker.io/library/python:3.12-slim
LABEL org.opencontainers.image.title="checkgm" \
      org.opencontainers.image.description="Score gut metagenomes against a healthy-adult reference baseline" \
      org.opencontainers.image.version="0.8.2" \
      org.opencontainers.image.source="https://github.com/khoa-yelo/Meta2"
ENV PIP_NO_CACHE_DIR=1 PYTHONUNBUFFERED=1 CHECKGM_BUNDLES=/bundles
WORKDIR /opt/checkgm
COPY pyproject.toml README.md LICENSE CHANGELOG.md ./
COPY src ./src
RUN pip install --upgrade pip && pip install . && rm -rf /root/.cache
VOLUME ["/bundles", "/work"]
WORKDIR /work
ENTRYPOINT ["checkgm"]
CMD ["--help"]

# refmb scoring image: the Python package only (no assembler or profiler). Use it to normalize and score MGnify v5 analyses,
# query directories or MetaPhlAn 3 / HUMAnN 3 tables against a bundle mounted into the container.
#   podman build -t refmb:0.8.1 -f container/Dockerfile .          # run from the refmb/ tool directory
#   (on a rootless host without subuid ranges add: --cgroup-manager=cgroupfs --events-backend=file)
#   podman run --rm -v /data/refmb_bundles:/bundles:ro -v $PWD:/work -e REFMB_BUNDLES=/bundles refmb:0.8.1 \
#       run-metaphlan --bundle gut-reads-adult-global-v0.2-lenient --input /work/S1.txt --reads-tsv /work/reads.tsv --out /work/S1_report
FROM docker.io/library/python:3.12-slim
LABEL org.opencontainers.image.title="refmb" \
      org.opencontainers.image.description="Score gut metagenomes against a healthy-adult reference baseline" \
      org.opencontainers.image.version="0.8.1" \
      org.opencontainers.image.source="https://github.com/khoa-yelo/Meta2"
ENV PIP_NO_CACHE_DIR=1 PYTHONUNBUFFERED=1 REFMB_BUNDLES=/bundles
WORKDIR /opt/refmb
COPY pyproject.toml README.md LICENSE CHANGELOG.md ./
COPY src ./src
COPY tests ./tests
RUN pip install --upgrade pip && pip install . && rm -rf /root/.cache
VOLUME ["/bundles", "/work"]
WORKDIR /work
ENTRYPOINT ["refmb"]
CMD ["--help"]

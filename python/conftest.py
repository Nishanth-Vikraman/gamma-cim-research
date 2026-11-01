"""Keep simulator test files inside the writable project workspace."""
import os
from pathlib import Path
import hashlib
import pytest


def pytest_configure(config):
    if config.option.basetemp is None:
        root=Path(__file__).resolve().parent/'results'/'.pytest'
        root.mkdir(parents=True,exist_ok=True)
        target=root/f'run_{os.getpid()}'
        assert target.parent == root
        config.option.basetemp=str(target)


@pytest.fixture
def tmp_path(tmp_path_factory,request):
    # LTspice XVII still encounters MAX_PATH limits for long test function names.
    name='t'+hashlib.sha256(request.node.nodeid.encode()).hexdigest()[:6]
    return tmp_path_factory.mktemp(name)

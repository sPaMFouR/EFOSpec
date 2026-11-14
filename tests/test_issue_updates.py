import json
from pathlib import Path
import pytest
from efospec.inventory import sha256
from efospec.review_handoff import issue_updates


@pytest.mark.parametrize('changed',[None,'raw','evidence','parent'])
def test_retry_notes_require_unchanged_parent_raw_and_evidence(tmp_path,changed):
    run=tmp_path/'reduction';run.mkdir();parent=run/'run.json';parent.write_text('{}')
    qc=tmp_path/'QC';qc.mkdir();raw=tmp_path/'raw.fits';raw.write_bytes(b'raw')
    evidence=qc/'retry.json';evidence.write_text(json.dumps(dict(parent_run=str(run),
        parent_sha256=sha256(parent),raw_sha256=sha256(raw),target='ATtest')))
    update=dict(object='ATtest',raw=str(raw),raw_sha256=sha256(raw),manifest=str(evidence),
        manifest_sha256=sha256(evidence),summary='Still unresolved after narrower bands.')
    (qc/'issue-updates.json').write_text(json.dumps(dict(parent_sha256=sha256(parent),updates=[update])))
    if changed:
        {'raw':raw,'evidence':evidence,'parent':parent}[changed].write_bytes(b'changed')
        with pytest.raises(ValueError):issue_updates(run)
    else:
        assert issue_updates(run)==[update]
        from efospec.epoch_report import EpochPage
        page=EpochPage(tmp_path)
        notice=page.failure_notice(run,dict(path=str(raw),object='ATtest',code='background_support_required',error='Old 30-pixel failure'))
        assert notice.index('Latest retry:')<notice.index('Original attempt')
        assert update['summary'] in notice and 'Old 30-pixel failure' in notice

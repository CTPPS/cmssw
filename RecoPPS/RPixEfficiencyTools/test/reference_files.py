"""
Files for reference efficiency analysis.
Taken from workflows/re-tracking-efficiency/python/reference_files.py
"""
import os
import re

def list_files(directory, pattern=None):
    """
    List files in directory. If pattern is provided (string or compiled re),
    only return files whose names match the regex (using re.search).
    """
    try:
        compiled = None
        if pattern is not None:
            compiled = pattern if hasattr(pattern, "search") else re.compile(pattern)
        files = []
        for f in os.listdir(directory):
            path = os.path.join(directory, f)
            if not os.path.isfile(path):
                continue
            if compiled is None or compiled.search(f):
                files.append(path)
        return files
    except OSError as e:
        print(f"Error accessing directory {directory}: {e}")
        return []

# List files in the specified directory
ref_files_2023 = list_files(
    '/eos/cms/store/group/dpg_ctpps/comm_ctpps/pps-automation/dev/pps-test-re-alcareco/zb_reference_2023_v2/366933/',
    pattern=r'^(?!.*plots).*\.root$'
)

ref_files_2024_preTS1 = list_files(
    '/eos/cms/store/group/dpg_ctpps/comm_ctpps/pps-automation/dev/pps-test-re-alcareco/2024_reference_preTS1_v1/378985/',
    pattern=r'^(?!.*plots).*\.root$'
)

ref_files_2024_postTS1 = list_files(
    '/eos/cms/store/group/dpg_ctpps/comm_ctpps/pps-automation/dev/pps-test-re-alcareco/2024_reference_postTS1_v1/382258/',
    pattern=r'^(?!.*plots).*\.root$'
)

REFERENCE_FILE_SETS = {
    "2023": ref_files_2023,
    "2024_preTS1": ref_files_2024_preTS1,
    "2024_postTS1": ref_files_2024_postTS1,
}
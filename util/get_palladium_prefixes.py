from get_id_to_paths import (
    get_uid_to_path__epic_and_hifi,
    get_prefixes
)

uid_to_path = get_uid_to_path__epic_and_hifi()
prefixes = get_prefixes(uid_to_path)

for prefix in prefixes:
    print(prefix)

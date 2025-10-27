from get_id_to_paths import (
    get_uid_to_path__epic_and_hifi,
    get_prefixes
)

def get_prefixes_wrapper():
    uid_to_path = get_uid_to_path__epic_and_hifi()
    prefixes = get_prefixes(uid_to_path)
    return prefixes

def main(): 
    prefixes = get_prefixes_wrapper()
    for prefix in prefixes:
        print(prefix)

if __name__ == '__main__': 
    main()
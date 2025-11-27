def prefix_columns(df, prefix, join_keys):
    # add prefixes to all columns except the join keys
    return df.rename({col: f"{prefix}_{col}" for col in df.columns if col not in join_keys})


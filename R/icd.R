# ========================
# ICD CODES - REGEX
# ========================

malig_insitu <- c(
    "^C(?!44)",           # Malignant neoplasms except skin cancer (C44)
    "^Z85(?!\\.828)",     # Personal history of malignant neoplasm except skin
    "^D0(?!4)[0-9]",      # Benign/in situ neoplasms except skin cancer (DO4)
    "^Z86.00",            # Personal history of in-situ neoplasm
    "^(?!173)(?:1[4-9][0-9]|2[0-3][0-9])",  # ICD-9 malignant/insitu neoplasms except 173
    "^V10(?!\\.83)"       # ICD-9 personal history of malignant neoplasm except skin
)

exclude_codes <- c("^C(?!44)", "^Z85(?!\\.828)", "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])", "^V10(?!\\.83)")

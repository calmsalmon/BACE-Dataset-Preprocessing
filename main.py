import pandas as pd
from rdkit import Chem

# Convert into df dataframe
df = pd.read_csv('bace.csv')

print(f"Original dataset size: {df.shape}")

# Remove empty columns
df = df.dropna(axis=1, how="all")

# Remove rows missing required data
required_columns = ["mol", "pIC50", "Class"]

df = df.dropna(subset=required_columns)

# Validate and canonicalize SMILES
def canonicalize_smiles(smiles):

    molecule = Chem.MolFromSmiles(str(smiles))

    # Invalid molecule
    if molecule is None:
        return None

    # Convert to standardized/canonical SMILES
    return Chem.MolToSmiles(
        molecule,
        canonical=True
    )


df["canonical_smiles"] = df["mol"].apply(
    canonicalize_smiles
)

# Remove invalid molecules
before_invalid_removal = len(df)

df = df.dropna(subset=["canonical_smiles"])

invalid_removed = (before_invalid_removal - len(df))

# Remove duplicate molecules
before_duplicate_removal = len(df)

df = df.drop_duplicates(
    subset=["canonical_smiles"],
    keep="first"
)

duplicates_removed = (before_duplicate_removal - len(df))

# Convert experiment data into numerical values
df["pIC50"] = pd.to_numeric(
    df["pIC50"],
    errors="coerce"
)

df["Class"] = pd.to_numeric(
    df["Class"],
    errors="coerce"
)

# Remove rows where experimental measurements
# could not be interpreted numerically.
df = df.dropna(
    subset=["pIC50", "Class"]
)

# These are the numerical molecular properties that
# can eventually be used by the classical portion
# of the hybrid algorithm
excluded_columns = {
    "CID",
    "pIC50",
    "Class"
}

descriptor_columns = [
    column
    for column in df.select_dtypes(
        include="number"
    ).columns
    if column not in excluded_columns
]

# A descriptor that has exactly the same value for
# every molecule provides no information.
constant_columns = [
    column
    for column in descriptor_columns
    if df[column].nunique(dropna=False) <= 1
]

df = df.drop(columns=constant_columns)

# Update descriptor list.
descriptor_columns = [
    column
    for column in df.select_dtypes(
        include="number"
    ).columns
    if column not in excluded_columns
]

# Use the median of each descriptor rather than
# deleting otherwise usable molecules.
for column in descriptor_columns:
    median_value = df[column].median()

    df[column] = df[column].fillna(
        median_value
    )

# Put the molecule and experimental ground truth
# at the beginning of the dataset.
important_columns = [
    "mol",
    "canonical_smiles",
    "pIC50",
    "Class"
]

remaining_columns = [
    column
    for column in df.columns
    if column not in important_columns
]

df = df[
    important_columns +
    remaining_columns
]

# Save data
output_file = "bace_preprocessed.csv"

df.to_csv(
    output_file,
    index=False
)

# Print summary
print("\n==============================")
print("PREPROCESSING COMPLETE")
print("==============================")

print(f"Final molecules: {len(df)}")
print(f"Final columns: {len(df.columns)}")

print(
    f"Invalid molecules removed: "
    f"{invalid_removed}"
)

print(
    f"Duplicate molecules removed: "
    f"{duplicates_removed}"
)

print(
    f"Descriptors retained: "
    f"{len(descriptor_columns)}"
)

print(
    f"\nPreprocessed dataset saved as: "
    f"{output_file}"
)

print("\nExperimental activity summary:")
print(
    df["pIC50"].describe()
)

print("\nClass distribution:")
print(
    df["Class"].value_counts()
)
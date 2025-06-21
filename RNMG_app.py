# Author: Kenneth Kusima
# Version 1.1
# 06/2025
# =============================================================================
# REACTION NETWORK MATRIX GENERATOR
# =============================================================================
# This app helps create stoichiometric and atomic matrices from chemical
# reaction mechanisms

# Import all the libraries we need
import streamlit as st
import pandas as pd
import numpy as np
import re
from io import StringIO
from collections import defaultdict, Counter
import math
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# =============================================================================
# PAGE SETUP AND CONFIGURATION
# =============================================================================
# Set up the page layout and basic configuration
st.set_page_config(
    page_title="RNMG - Reaction Network Matrix Generator",
    page_icon="🔗",
    layout="wide"  # Uses full width of browser
)

# Main title and description
st.title("RNMG - Reaction Network Matrix Generator")
st.markdown("Create stoichiometric and atomic matrices from chemical reaction mechanisms with equilibrium analysis")

# =============================================================================
# SESSION STATE INITIALIZATION
# =============================================================================
# Streamlit uses session state to remember data between runs
# We initialize all our variables here so they persist
if 'reactions' not in st.session_state:
    st.session_state.reactions = []  # Stores all the reactions we've added
if 'species_list' not in st.session_state:
    st.session_state.species_list = []  # List of all unique species
if 'atomic_matrix' not in st.session_state:
    st.session_state.atomic_matrix = pd.DataFrame()  # A matrix showing atomic composition
if 'stoich_matrix' not in st.session_state:
    st.session_state.stoich_matrix = pd.DataFrame()  # Matrix showing reaction stoichiometry
if 'param_data' not in st.session_state:
    st.session_state.param_data = []  # Parameters like rate constants and pressures
if 'matrix_consistency' not in st.session_state:
    st.session_state.matrix_consistency = {'consistent': True, 'atomic_species': [], 'stoich_species': []}
if 'mass_balance_result' not in st.session_state:
    st.session_state.mass_balance_result = {'balanced': None, 'message': ''}
if 'app_mode' not in st.session_state:
    st.session_state.app_mode = 'Welcome'  # Can be 'Welcome', 'Generator' or 'Upload'
if 'uploaded_files' not in st.session_state:
    st.session_state.uploaded_files = {'atomic': None, 'stoich': None, 'params': None}
if 'analysis_mode' not in st.session_state:
    st.session_state.analysis_mode = 'Mass Balance'  # Type of analysis
if 'sigma_values' not in st.session_state:
    st.session_state.sigma_values = []  # Stoichiometric numbers for thermodynamics
if 'equilibrium_constants' not in st.session_state:
    st.session_state.equilibrium_constants = []  # Calculated equilibrium constants
if 'equilibrium_results' not in st.session_state:
    st.session_state.equilibrium_results = {}  # Equilibrium analysis results

# =============================================================================
# HELPER FUNCTIONS FOR PARSING REACTIONS AND SPECIES
# =============================================================================

def parse_species(species_str):
    """
    Takes a species string like "2O*" or "CO(g)" and figures out
    what the actual species is and how many of them we have.
    Returns the species name and its coefficient.
    """
    species_str = species_str.strip()  # Remove any extra spaces
    # Look for numbers at the beginning (like "2" in "2O*")
    coeff_match = re.match(r'^(\d+)(.+)$', species_str)
    if coeff_match:
        coeff = int(coeff_match.group(1))  # The number part
        species = coeff_match.group(2)  # The species part
    else:
        coeff = 1  # If no number, coefficient is 1
        species = species_str
    return species, coeff

def parse_reaction(reaction_str):
    """
    Takes a reaction string like "CO(g) + * ⇌ CO*" and breaks it down
    into reactants and products with their coefficients.
    This is the heart of the reaction parsing system!
    """
    # We only accept the equilibrium arrow ⇌ to keep things consistent
    if '⇌' not in reaction_str:
        st.error(f"❌ Please use only ⇌ arrows. Copy the arrow from the box above. Your input: {reaction_str}")
        return None, None, None

    # Split the reaction at the arrow
    parts = reaction_str.split('⇌')
    if len(parts) != 2:
        st.error(f"Invalid reaction format: {reaction_str}")
        return None, None, None

    reactants_str, products_str = parts

    # Parse the reactants side
    reactants = {}
    for reactant in reactants_str.split('+'):
        reactant = reactant.strip()
        if reactant:  # Skip empty strings
            species, coeff = parse_species(reactant)
            # If we see the same species twice, add up the coefficients
            reactants[species] = reactants.get(species, 0) + coeff

    # Parse the products side (same logic as reactants)
    products = {}
    for product in products_str.split('+'):
        product = product.strip()
        if product:
            species, coeff = parse_species(product)
            products[species] = products.get(species, 0) + coeff

    return reactants, products, reaction_str

# =============================================================================
# CHEMICAL FORMULA PARSING FUNCTIONS
# =============================================================================

def get_atomic_composition():
    """
    Provides a fallback dictionary for surface site types.
    These are the basic surface sites that can't be parsed automatically.
    """
    compositions = {
        '*': {'*': 1},  # Standard surface site
        '_': {'_': 1},  # Alternative surface site type
        '#': {'#': 1},  # Alternative surface site type
        '@': {'@': 1},  # Alternative surface site type
        '&': {'&': 1},  # Alternative surface site type
    }
    return compositions

def parse_species_formula(formula):
    """
    This is our comprehensive chemical formula parser! It can handle complex
    formulas like Ca(OH)2, Al2(SO4)3, and surface species with various site types.
    Much more powerful than the old hardcoded approach.
    """
    composition = {}
    # Define all possible surface site symbols
    surface_symbols = ['*', '_', '#', '@', '&']

    # Handle surface site cases
    for symbol in surface_symbols:
        if formula == symbol:
            return {symbol: 1}

    # Check if this is an adsorbed species (ends with any surface symbol)
    surface_site = None
    for symbol in surface_symbols:
        if formula.endswith(symbol):
            surface_site = symbol
            formula = formula[:-1]  # Remove the surface symbol for now
            break

    def expand_parentheses(formula_str):
        """
        Handles parentheses in chemical formulas like Ca(OH)2.
        This expands Ca(OH)2 into CaO2H2 so we can parse it easily.
        """
        import re
        # Keep expanding until no more parentheses
        while '(' in formula_str:
            # Find the innermost parentheses (handles nested cases)
            match = re.search(r'\(([^()]+)\)(\d*)', formula_str)
            if not match:
                break  # Something went wrong, bail out

            group_content = match.group(1)  # What's inside the parentheses
            multiplier = int(match.group(2)) if match.group(2) else 1  # Number after parentheses

            # Expand everything inside the parentheses
            expanded = ''
            group_elements = re.findall(r'([A-Z][a-z]?)(\d*)', group_content)
            for element, count in group_elements:
                count = int(count) if count else 1
                expanded += element + str(count * multiplier)

            # Replace the parentheses group with the expanded version
            formula_str = formula_str[:match.start()] + expanded + formula_str[match.end():]

        return formula_str

    # First, expand any parentheses we find
    expanded_formula = expand_parentheses(formula)

    # Now parse all the element-number pairs
    import re
    pattern = r'([A-Z][a-z]?)(\d*)'  # Matches things like "Ca", "O2", "Al", etc.
    matches = re.findall(pattern, expanded_formula)

    # Build up our composition dictionary
    for element, count in matches:
        count = int(count) if count else 1  # If no number, assume 1
        composition[element] = composition.get(element, 0) + count

    # Add the surface site back if this was an adsorbed species
    if surface_site:
        composition[surface_site] = 1

    # If we couldn't parse anything, mark it as unknown
    if not composition and formula:
        composition = {'Unknown': 1}
        print(f"Warning: Could not parse formula '{formula}'. Treating as unknown species.")

    return composition

def get_species_atomic_composition(species):
    """
    Main function that gets the atomic composition for any species.
    First tries the fallback dictionary, then uses our smart parser.
    """
    # Remove phase notation like (g), (l), (s) - we don't need it for composition
    clean_species = re.sub(r'\([gls]\)', '', species)

    # Check the basic fallback dictionary first
    compositions = get_atomic_composition()
    if clean_species in compositions:
        return compositions[clean_species]

    # Use our comprehensive parser for everything else
    return parse_species_formula(clean_species)

# =============================================================================
# PARAMETER AND MATRIX GENERATION FUNCTIONS
# =============================================================================

def create_default_parameters():
    """
    Creates a default set of parameters for the reactions we've added.
    This includes temperature, gas constant, pressures, rate constants,
    and some placeholder constants. Users can edit these later.
    """
    params = []

    # Global thermodynamic parameters
    params.append({'Reaction_Descrp': '', 'Parameter': 'T', 'Values': 320.0, 'Units': 'K'})
    params.append({'Reaction_Descrp': '', 'Parameter': 'R', 'Values': 8.31446, 'Units': 'JK^-1mol^-1'})

    # Create pressure parameters for all gas phase species
    surface_symbols = ['*', '_', '#', '@', '&']
    gas_species = [s for s in st.session_state.species_list if not any(s.endswith(sym) for sym in surface_symbols)]

    for i, species in enumerate(gas_species):
        params.append({'Reaction_Descrp': species, 'Parameter': f'P{i+1}', 'Values': 1.0e-8, 'Units': 'bar'})

    # Create forward and reverse rate constants for each reaction
    for i in range(len(st.session_state.reactions)):
        reaction_id = f'r{i+1}'
        params.append({'Reaction_Descrp': reaction_id, 'Parameter': f'k{i+1}f', 'Values': 1.0, 'Units': '1/s'})
        params.append({'Reaction_Descrp': '', 'Parameter': f'k{i+1}r', 'Values': 1.0, 'Units': '1/s'})

    # Add some placeholder constants (users might need these for coverage dependence in kinetic modeling)
    for i in range(len(st.session_state.reactions)):
        for direction in ['f', 'r']:  # Forward and reverse
            for const_type in ['a', 'b', 'c']:  # Different types of constants
                params.append({'Reaction_Descrp': 'const', 'Parameter': f'{const_type}{i+1}{direction}', 'Values': 1.0, 'Units': '-'})

    return params

def update_matrices():
    """
    This is a big function that rebuilds both the atomic and stoichiometric
    matrices whenever we add or remove reactions. It makes sure everything
    stays consistent and properly ordered.
    """
    # If no reactions, just clear everything
    if not st.session_state.reactions:
        st.session_state.atomic_matrix = pd.DataFrame()
        st.session_state.stoich_matrix = pd.DataFrame()
        st.session_state.species_list = []
        return

    # Collect all unique species from all reactions
    all_species = set()
    for _, reactants, products, _ in st.session_state.reactions:
        all_species.update(reactants.keys())
        all_species.update(products.keys())

    # Define surface symbols for categorization
    surface_symbols = ['*', '_', '#', '@', '&']

    # Organize species: gas phase first, then surface species, with bare sites at the end
    gas_species = sorted([s for s in all_species if not any(s.endswith(sym) for sym in surface_symbols)])
    surface_species = [s for s in all_species if any(s.endswith(sym) for sym in surface_symbols)]

    # Sort surface species but keep bare sites at the end
    bare_sites = [s for s in surface_species if s in surface_symbols]
    other_surface = sorted([s for s in surface_species if s not in surface_symbols])
    surface_species = other_surface + bare_sites

    # Create our final ordered species list
    ordered_species = gas_species + surface_species
    st.session_state.species_list = ordered_species

    # ======================
    # BUILD THE ATOMIC MATRIX
    # ======================
    # Find all the different atoms we need to track
    all_atoms = set()
    for species in ordered_species:
        composition = get_species_atomic_composition(species)
        all_atoms.update(composition.keys())

    # Sort atoms with surface sites at the end
    other_atoms = sorted([atom for atom in all_atoms if atom not in surface_symbols])
    surface_atoms = [atom for atom in all_atoms if atom in surface_symbols]
    all_atoms = other_atoms + surface_atoms

    # Build the atomic matrix data
    atomic_data = []
    for atom in all_atoms:
        row = [atom]  # First column is the atom name
        for species in ordered_species:
            composition = get_species_atomic_composition(species)
            row.append(composition.get(atom, 0))  # How many of this atom in this species
        atomic_data.append(row)

    # Create the atomic matrix DataFrame
    columns = ['A\\S'] + ordered_species  # A\S means "Atoms \ Species"
    st.session_state.atomic_matrix = pd.DataFrame(atomic_data, columns=columns)

    # ============================
    # BUILD THE STOICHIOMETRIC MATRIX
    # ============================
    stoich_data = []

    # Create column headers with proper notation
    gas_species = [s for s in ordered_species if not any(s.endswith(sym) for sym in surface_symbols)]
    surface_species = [s for s in ordered_species if any(s.endswith(sym) for sym in surface_symbols)]

    stoich_columns = ['r\\S'] + [f'P_{s}' for s in gas_species] + [f'theta_{s}' for s in surface_species]

    # Build each row (one for each reaction)
    for i, (_, reactants, products, _) in enumerate(st.session_state.reactions):
        row = [f'r{i+1}']  # First column is reaction name

        # Gas phase species columns
        for species in gas_species:
            # Net stoichiometric coefficient = products - reactants
            net_coeff = products.get(species, 0) - reactants.get(species, 0)
            row.append(net_coeff)

        # Surface species columns
        for species in surface_species:
            net_coeff = products.get(species, 0) - reactants.get(species, 0)
            row.append(net_coeff)

        stoich_data.append(row)

    # Create the stoichiometric matrix DataFrame
    st.session_state.stoich_matrix = pd.DataFrame(stoich_data, columns=stoich_columns)

    # Check that both matrices have consistent species ordering
    verify_matrix_consistency()

def verify_matrix_consistency():
    """
    Double-checks that our atomic and stoichiometric matrices have
    the same species in the same order. This is crucial for mass balance
    calculations to work properly.
    """
    if st.session_state.atomic_matrix.empty or st.session_state.stoich_matrix.empty:
        return True

    # Get species from atomic matrix (skip first column which is atom names)
    atomic_species = st.session_state.atomic_matrix.columns[1:].tolist()

    # Get species from stoichiometric matrix and clean up the prefixes
    stoich_species = []
    for col in st.session_state.stoich_matrix.columns[1:]:
        if col.startswith('P_'):
            stoich_species.append(col[2:])  # Remove 'P_' prefix
        elif col.startswith('theta_'):
            stoich_species.append(col[6:])  # Remove 'theta_' prefix

    # Check if the lists match
    consistency_check = atomic_species == stoich_species

    # Store the result for display purposes
    st.session_state.matrix_consistency = {
        'consistent': consistency_check,
        'atomic_species': atomic_species,
        'stoich_species': stoich_species
    }

    return consistency_check

# =============================================================================
# ANALYSIS AND CALCULATION FUNCTIONS
# =============================================================================

def check_mass_balance(atomic_matrix=None, stoich_matrix=None):
    """
    Performs mass balance checking using the fundamental equation A*v=0,
    where A is the atomic matrix and v is the stoichiometric vector.
    This checks if atoms are conserved across all reactions.
    """
    # Use provided matrices or default to session state
    if atomic_matrix is None:
        atomic_matrix = st.session_state.atomic_matrix
    if stoich_matrix is None:
        stoich_matrix = st.session_state.stoich_matrix

    if atomic_matrix.empty or stoich_matrix.empty:
        return False, "Matrices not available for mass balance check"

    # Extract the numerical parts (skip label columns)
    atomic_matrix_vals = atomic_matrix.iloc[:, 1:].values
    stoich_matrix_vals = stoich_matrix.iloc[:, 1:].values

    errors = []
    mass_balanced = True

    # Check each reaction individually
    for i, reaction_row in enumerate(stoich_matrix_vals):
        # Calculate A * v for this reaction
        mass_balance_result = np.dot(atomic_matrix_vals, reaction_row)

        # If any element is not close to zero, we have an imbalance
        if not np.allclose(mass_balance_result, 0, atol=1e-10):
            mass_balanced = False
            reaction_name = stoich_matrix.iloc[i, 0]
            errors.append({
                'reaction': i + 1,
                'reaction_name': reaction_name,
                'imbalance': mass_balance_result.tolist()
            })

    # Prepare the result message
    if mass_balanced:
        return True, "✅ Mass is conserved across all reactions."
    else:
        error_message = "❌ Mass is NOT conserved in the following reactions:\n"
        for error in errors:
            error_message += f"• Reaction {error['reaction']} ({error['reaction_name']}): Imbalance = {error['imbalance']}\n"
        error_message += "\nPlease check and correct the Atomic and/or Stoichiometric matrices."
        return False, error_message

def calculate_equilibrium_constants(params_df, sigma_values, k_pattern="auto"):
    """
    Comprehensive function to calculate equilibrium constants from rate constants with comprehensive error handling.
    This function calculates individual K_eq values and the overall equilibrium constant using K_overall = ∏K_i^σ_i

    Args:
        params_df: DataFrame with parameters including rate constants
        sigma_values: List of stoichiometric numbers for each reaction
        k_pattern: How rate constants are organized ("auto" or "alternating")

    Returns:
        Tuple of (K_eq_list, K_overall, temperature, gas_constant, k_forward, k_reverse, ln_K_list, gibbs_free_energy_list)
    """
    try:
        # Extract temperature and gas constant from parameters with better error handling
        T_row = params_df[params_df['Parameter'].str.lower() == 't']
        R_row = params_df[params_df['Parameter'].str.lower() == 'r']

        # Use defaults if not found, with user warnings
        if T_row.empty:
            st.warning("Temperature (T) not found in parameters, using default 320 K")
            T = 320.0
        else:
            T = float(T_row['Values'].iloc[0])

        if R_row.empty:
            st.warning("Gas constant (R) not found in parameters, using default 8.314 J/(mol·K)")
            R = 8.314
        else:
            R = float(R_row['Values'].iloc[0])

        # Find all rate constants with improved pattern matching
        k_forward = []
        k_reverse = []

        # Get all parameters that look like rate constants
        k_params = params_df[params_df['Parameter'].str.contains(r'k', regex=True, na=False, case=False)]
        k_params = k_params.sort_values('Parameter')

        if k_pattern == "alternating":
            # Every other k is reverse (k1=forward, k2=reverse, etc.)
            k_values = []
            for _, row in k_params.iterrows():
                param_name = str(row['Parameter'])
                if 'k' in param_name.lower():
                    try:
                        k_values.append(float(row['Values']))
                    except (ValueError, TypeError):
                        st.warning(f"Could not parse rate constant {param_name}, using default 1.0")
                        k_values.append(1.0)

            for i in range(0, len(k_values), 2):
                if i < len(k_values):
                    k_forward.append(k_values[i])
                    if i + 1 < len(k_values):
                        k_reverse.append(k_values[i + 1])
                    else:
                        k_reverse.append(1.0)  # Default if missing
        else:
            # Try to automatically detect the naming pattern
            for i in range(len(sigma_values)):
                # Try different possible naming conventions
                possible_forward = [f'k{i+1}f', f'k{i+1}F', f'k{i+1}_f', f'k{i+1}', f'k{2*i+1}']
                possible_reverse = [f'k{i+1}r', f'k{i+1}R', f'k{i+1}_r', f'k{i+1}_rev', f'k{2*i+2}']

                kf_found = False
                kr_found = False

                # Look for forward rate constant
                for kf_name in possible_forward:
                    kf_row = params_df[params_df['Parameter'] == kf_name]
                    if not kf_row.empty:
                        try:
                            k_forward.append(float(kf_row['Values'].iloc[0]))
                            kf_found = True
                            break
                        except (ValueError, TypeError):
                            continue

                # Look for reverse rate constant
                for kr_name in possible_reverse:
                    kr_row = params_df[params_df['Parameter'] == kr_name]
                    if not kr_row.empty:
                        try:
                            k_reverse.append(float(kr_row['Values'].iloc[0]))
                            kr_found = True
                            break
                        except (ValueError, TypeError):
                            continue

                # Use defaults if we couldn't find them
                if not kf_found:
                    st.warning(f"Forward rate constant for reaction {i+1} not found, using default 1.0")
                    k_forward.append(1.0)
                if not kr_found:
                    st.warning(f"Reverse rate constant for reaction {i+1} not found, using default 1.0")
                    k_reverse.append(1.0)

        # Make sure we have the right number of rate constants
        while len(k_forward) < len(sigma_values):
            k_forward.append(1.0)
        while len(k_reverse) < len(sigma_values):
            k_reverse.append(1.0)

        # Calculate individual equilibrium constants with error checking
        K_eq = []
        ln_K_list = []
        gibbs_free_energy_list = []

        for i in range(len(sigma_values)):
            if k_reverse[i] != 0 and not np.isnan(k_reverse[i]) and not np.isinf(k_reverse[i]):
                K_i = k_forward[i] / k_reverse[i]
                if np.isnan(K_i) or np.isinf(K_i):
                    st.warning(f"Invalid equilibrium constant for reaction {i+1}, using 1.0")
                    K_i = 1.0
                K_eq.append(K_i)

                # Calculate ln(K)
                if K_i > 0:
                    ln_K_i = np.log(K_i)
                    ln_K_list.append(ln_K_i)
                    # Calculate Gibbs free energy: ΔG = -RT ln(K)
                    delta_G = -R * T * ln_K_i
                    gibbs_free_energy_list.append(delta_G)
                else:
                    ln_K_list.append(0.0)
                    gibbs_free_energy_list.append(0.0)
            else:
                st.warning(f"Zero or invalid reverse rate constant for reaction {i+1}, using K_eq = 1.0")
                K_eq.append(1.0)
                ln_K_list.append(0.0)
                gibbs_free_energy_list.append(0.0)

        # Calculate overall equilibrium constant: K_overall = ∏(K_i^σ_i)
        K_overall = 1.0
        for i, (K_i, sigma_i) in enumerate(zip(K_eq, sigma_values)):
            if K_i > 0 and not np.isnan(K_i) and not np.isinf(K_i):  # Make sure we can take the power
                try:
                    contribution = K_i ** sigma_i
                    if np.isnan(contribution) or np.isinf(contribution):
                        st.warning(f"Invalid contribution from reaction {i+1}, skipping")
                        continue
                    K_overall *= contribution
                except (OverflowError, ZeroDivisionError):
                    st.warning(f"Mathematical error in reaction {i+1}, using contribution = 1.0")
                    continue
            else:
                st.warning(f"Invalid equilibrium constant for reaction {i+1} in overall calculation")

        # Final validation of K_overall
        if np.isnan(K_overall) or np.isinf(K_overall):
            st.error("Overall equilibrium constant is invalid, using 1.0")
            K_overall = 1.0

        # Calculate overall ln(K) and Gibbs free energy
        if K_overall > 0:
            ln_K_overall = np.log(K_overall)
            delta_G_overall = -R * T * ln_K_overall
        else:
            ln_K_overall = 0.0
            delta_G_overall = 0.0

        return K_eq, K_overall, T, R, k_forward, k_reverse, ln_K_list, gibbs_free_energy_list, ln_K_overall, delta_G_overall

    except Exception as e:
        st.error(f"Error in equilibrium constant calculation: {str(e)}")
        # Return safe defaults if calculation fails
        n_reactions = len(sigma_values)
        return [1.0] * n_reactions, 1.0, 320.0, 8.314, [1.0] * n_reactions, [1.0] * n_reactions, [0.0] * n_reactions, [0.0] * n_reactions, 0.0, 0.0

def create_equilibrium_plots(K_eq_list, sigma_values, reaction_names, K_overall):
    """
    Creates comprehensive plots for equilibrium analysis
    """
    # Individual equilibrium constants bar chart
    fig_individual = go.Figure()

    # Color code based on equilibrium position
    colors = []
    for K_i in K_eq_list:
        if K_i > 10:
            colors.append('darkgreen')  # Strongly forward
        elif K_i > 1:
            colors.append('lightgreen')  # Forward
        elif K_i == 1:
            colors.append('yellow')  # Balanced
        elif K_i > 0.1:
            colors.append('orange')  # Reverse
        else:
            colors.append('red')  # Strongly reverse

    fig_individual.add_trace(go.Bar(
        x=reaction_names,
        y=K_eq_list,
        marker_color=colors,
        hovertemplate="%{x}<br>K_eq: %{y:.2e}<br>σ: %{customdata:.3f}",
        customdata=sigma_values,
        name="Individual K_eq"
    ))

    fig_individual.update_layout(
        title="Individual Equilibrium Constants (K_eq = k_forward / k_reverse)",
        xaxis_title="Reaction",
        yaxis_title="K_eq",
        yaxis_type="log",
        showlegend=False,
        height=500
    )

    # Contribution to overall equilibrium
    contributions = [K_i ** sigma_i for K_i, sigma_i in zip(K_eq_list, sigma_values)]

    fig_contribution = go.Figure()

    fig_contribution.add_trace(go.Bar(
        x=reaction_names,
        y=contributions,
        marker_color='steelblue',
        hovertemplate="%{x}<br>K_eq^σ: %{y:.2e}<br>K_eq: %{customdata[0]:.2e}<br>σ: %{customdata[1]:.3f}",
        customdata=list(zip(K_eq_list, sigma_values)),
        name="K_eq^σ"
    ))

    fig_contribution.update_layout(
        title=f"Contribution to Overall Equilibrium (K_overall = {K_overall:.2e})",
        xaxis_title="Reaction",
        yaxis_title="K_eq^σ",
        yaxis_type="log",
        showlegend=False,
        height=500
    )

    return fig_individual, fig_contribution

def perform_analysis(atomic_matrix, stoich_matrix, params_df):
    """
    Unified analysis function that can be called from both modes
    """
    st.subheader("🔍 Analysis Options")

    # Analysis type selection
    analysis_type = st.selectbox(
        "Choose the type of analysis to perform:",
        ["Mass Balance",
         "Equilibrium Analysis"]
    )

    # Extract the actual analysis type
    if "Mass Balance" in analysis_type:
        selected_analysis = "Mass Balance"
    else:
        selected_analysis = "Equilibrium Analysis"

    if selected_analysis == "Mass Balance":
        st.subheader("⚖️ Mass Balance Analysis")
        st.caption("Checks if atoms are conserved across all reactions using A × ν = 0")

        if not atomic_matrix.empty and not stoich_matrix.empty:
            if st.button("Check Mass Balance", type="primary"):
                balanced, message = check_mass_balance(atomic_matrix, stoich_matrix)
                if balanced:
                    st.success(message)
                else:
                    st.error(message)
        else:
            st.warning("Please provide both atomic and stoichiometric matrices to perform mass balance analysis")

    elif selected_analysis == "Equilibrium Analysis":
        st.subheader("🔬 Equilibrium Analysis")
        st.caption("Calculate equilibrium constants from kinetic parameters")

        if not stoich_matrix.empty and not params_df.empty:
            # Input for stoichiometric numbers
            st.write("**Stoichiometric Numbers (σ) for each reaction:**")
            st.caption("These determine how each reaction contributes to the overall equilibrium constant")

            n_reactions = len(stoich_matrix) if not stoich_matrix.empty else 0

            if n_reactions > 0:
                sigma_input = st.text_input(
                    f"Enter {n_reactions} stoichiometric numbers (comma-separated):",
                    value=",".join(["1"] * n_reactions),
                    help="Example: 1,1,-1,1 for 4 reactions. Use negative values to reverse reaction contribution."
                )

                try:
                    sigma_values = [float(x.strip()) for x in sigma_input.split(',')]
                    if len(sigma_values) != n_reactions:
                        st.error(f"Please enter exactly {n_reactions} values")
                    else:
                        st.session_state.sigma_values = sigma_values

                        if st.button("Calculate Equilibrium Constants", type="primary"):
                            K_eq_list, K_overall, T, R, k_forward, k_reverse, ln_K_list, gibbs_free_energy_list, ln_K_overall, delta_G_overall = calculate_equilibrium_constants(
                                params_df, sigma_values
                            )

                            # Store results
                            st.session_state.equilibrium_results = {
                                'K_eq_list': K_eq_list,
                                'K_overall': K_overall,
                                'T': T,
                                'R': R,
                                'k_forward': k_forward,
                                'k_reverse': k_reverse,
                                'sigma_values': sigma_values,
                                'ln_K_list': ln_K_list,
                                'gibbs_free_energy_list': gibbs_free_energy_list,
                                'ln_K_overall': ln_K_overall,
                                'delta_G_overall': delta_G_overall
                            }

                            # Display overall results
                            st.success(f"**Overall Equilibrium Constant:** K_overall = {K_overall:.2e}")
                            st.info(f"**Overall ln(K):** ln(K_overall) = {ln_K_overall:.3f}")
                            st.info(f"**Gas Constant:** R = 8.314 J/mol K")
                            st.info(f"**Overall Gibbs Free Energy:** ΔG_overall = -RTlnK = {delta_G_overall*1e-3:.2f} kJ/mol")

                            # Results table
                            results_data = []
                            reaction_names = [f"r{i+1}" for i in range(len(K_eq_list))]

                            for i, (K_eq, sigma, kf, kr, ln_K, delta_G) in enumerate(zip(K_eq_list, sigma_values, k_forward, k_reverse, ln_K_list, gibbs_free_energy_list)):
                                results_data.append({
                                    'Reaction': reaction_names[i],
                                    'k_forward': f"{kf:.2e}",
                                    'k_reverse': f"{kr:.2e}",
                                    'K_eq': f"{K_eq:.2e}",
                                    'ln(K)': f"{ln_K:.3f}",
                                    'ΔG (J/mol)': f"{delta_G:.2f}",
                                    'σ': f"{sigma:.3f}",
                                    'K_eq^σ': f"{K_eq**sigma:.2e}"
                                })

                            results_df = pd.DataFrame(results_data)
                            st.dataframe(results_df, use_container_width=True)

                            # Create and display plots
                            fig_individual, fig_contribution = create_equilibrium_plots(
                                K_eq_list, sigma_values, reaction_names, K_overall
                            )

                            st.plotly_chart(fig_individual, use_container_width=True)
                            st.plotly_chart(fig_contribution, use_container_width=True)

                            # Interpretation
                            st.subheader("📖 Interpretation")
                            st.markdown("""
                            **Natural Logarithm of Equilibrium Constant (ln K):**
                            - **ln K > 0**: Reaction favors products
                            - **ln K = 0**: Reaction is balanced
                            - **ln K < 0**: Reaction favors reactants

                            **Gibbs Free Energy (ΔG):**
                            - **ΔG < 0**: Reaction is thermodynamically favorable (spontaneous)
                            - **ΔG = 0**: Reaction is at equilibrium
                            - **ΔG > 0**: Reaction is thermodynamically unfavorable (non-spontaneous)

                            **Overall Equilibrium Constant:**
                            K_overall = ∏(K_i^σ_i) represents the combined thermodynamic effect of all reactions
                            """)

                except ValueError:
                    st.error("Please enter valid numbers separated by commas")
        else:
            st.warning("Please provide stoichiometric matrix and parameters to perform equilibrium analysis")

# =============================================================================
# SIDEBAR CONFIGURATION
# =============================================================================

# Sidebar for mode selection
st.sidebar.title("🔧 Mode Selection")
st.sidebar.markdown("Choose your working mode:")

app_mode = st.sidebar.selectbox(
    "Select Mode:",
    ["🏠 Welcome", "🔬 Generator Mode", "📁 Upload Mode"],
    index=0 if st.session_state.app_mode == 'Welcome' else (1 if st.session_state.app_mode == 'Generator' else 2),
    help="Welcome: Introduction | Generator: Create matrices from reactions | Upload: Analyze existing files"
)

# Extract the actual mode from the selection
if "Welcome" in app_mode:
    selected_mode = "Welcome"
elif "Generator" in app_mode:
    selected_mode = "Generator"
else:
    selected_mode = "Upload"

# Update session state
st.session_state.app_mode = selected_mode

# Add some information in the sidebar
st.sidebar.markdown("---")
st.sidebar.markdown("### 📋 Current Mode")
if selected_mode == "Welcome":
    st.sidebar.info("🏠 **Welcome Screen**")
    st.sidebar.markdown("Learn about RNMG capabilities")
elif selected_mode == "Generator":
    st.sidebar.success("🔬 **Generator Mode**")
    st.sidebar.markdown("Build matrices from chemical reactions")
else:
    st.sidebar.info("📁 **Upload Mode**")
    st.sidebar.markdown("Analyze existing matrix files")

# Add help section in sidebar
st.sidebar.markdown("---")
st.sidebar.markdown("### ❓ Quick Help")
if selected_mode == "Welcome":
    st.sidebar.markdown("""
    **Welcome to RNMG!**
    
    Choose a mode to get started:
    - **Generator**: Build from reactions
    - **Upload**: Analyze existing files
    """)
elif selected_mode == "Generator":
    st.sidebar.markdown("""
    **Steps:**
    1. Enter reactions using ⇌ arrow
    2. Click 'Add Reaction'
    3. Generate parameters
    4. Perform analysis
    """)
else:
    st.sidebar.markdown("""
    **Steps:**
    1. Upload CSV files
    2. Review matrices
    3. Perform analysis
    4. Download results
    """)

# =============================================================================
# MAIN APPLICATION INTERFACE
# =============================================================================

st.markdown("---")

# Welcome Screen
if selected_mode == "Welcome":
    # Welcome content with comprehensive introduction
    st.markdown("""
    ## 🎉 Welcome to RNMG!
    
    **RNMG (Reaction Network Matrix Generator)** is a tool for researchers working with reaction mechanisms. This application helps you create and analyze stoichiometric and atomic matrices from chemical reactions.
    """)
    
    # Feature overview
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        ### 🔬 **Generator Mode**
        
        **Build reaction networks:**
        
        ✅ **Input reactions one by one**  
        ✅ **Automatic matrix generation**  
        ✅ **Smart species ordering**  
        ✅ **Parameter file creation**  
        ✅ **Real-time validation**  
        
        **Supported Features:**
        - Gas-solid heterogeneous catalysis
        - Multiple surface site types (`*`, `_`, `#`, `@`, `&`)
        - Complex chemical formulas (Ca(OH)₂, Al₂(SO₄)₃)
        - Automatic mass balance checking
        """)
    
    with col2:
        st.markdown("""
        ### 📁 **Upload Mode**
        
        **Analyzing pre-made matrices:**
        
        ✅ **Upload CSV matrix files**  
        ✅ **Import parameter sets**  
        ✅ **Equilibrium Analysis**  
        ✅ **Thermodynamic calculations**  
        ✅ **Visual data exploration**  
        
        **Analysis Capabilities:**
        - Mass balance verification (A × ν = 0)
        - Equilibrium constant calculations
        - Gibbs free energy analysis
        - Interactive plotting and visualization
        """)
    
    # What RNMG generates
    st.markdown("""
    ### 📊 **What RNMG Generates**
    
    **Three essential files for reaction network analysis:**
    """)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("""
        **🧮 Atomic Matrix (A)**
        - Shows atomic composition of each species
        - Rows: Elements (C, H, O, *, etc.)
        - Columns: Chemical species
        - Used for mass balance verification
        """)
    
    with col2:
        st.markdown("""
        **⚖️ Stoichiometric Matrix (ν)**
        - Net coefficients for each reaction
        - Rows: Reactions (r1, r2, r3, etc.)
        - Columns: Species (gas and surface)
        - Foundation for kinetic modeling
        """)
    
    with col3:
        st.markdown("""
        **⚙️ Parameter File**
        - Temperature and gas constant
        - Pressure values for gas species
        - Forward/reverse rate constants
        - Customizable constants for coverage dependence (Bragg Williams)
        """)
    
    # Analysis capabilities
    st.markdown("""
    ### 🔍 **Analysis Capabilities**
    
    **RNMG provides powerful analysis tools for your reaction networks:**
    """)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        **⚖️ Mass Balance Analysis**
        - Verifies atom conservation: **A × ν = 0**
        - Identifies problematic reactions
        - Essential for mechanism validation
        - Automatic error detection and reporting
        """)
    
    with col2:
        st.markdown("""
        **🔬 Equilibrium Analysis**
        - Calculate individual equilibrium constants: **K = k_f/k_r**
        - Overall equilibrium: **K_overall = ∏K_i^σ_i**
        - Gibbs free energy calculations: **ΔG = -RT ln(K)**
        - Interactive plots and thermodynamic insights
        """)
    
    # Important rules and limitations
    st.markdown("""
    ### ⚠️ **Important Rules & Guidelines**
    
    **For best results, please follow these guidelines:**
    """)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        **✅ DO:**
        - Use **only ⇌ arrows** (copy from input box)
        - Enter **ONE reaction at a time**
        - Use supported surface sites: `*`, `_`, `#`, `@`, `&`
        - Check mass balance before analysis
        - Verify matrix consistency
        """)
    
    with col2:
        st.markdown("""
        **❌ DON'T:**
        - Mix different arrow types (→, ←, ↔)
        - Enter multiple reactions in one line
        - Use unsupported surface site symbols
        - Skip mass balance verification
        - Ignore consistency warnings
        """)
    
    # Current limitations
    st.markdown("""
    ### 📝 **Current Limitations**
    
    **RNMG is currently optimized for:**
    - **Gas-solid heterogeneous catalysis** (gas-surface reactions)
    - **Single and Two phase surface reactions**
    - **Elementary reaction steps**
    """)
    
    # Getting started
    st.markdown("""
    ### 🚀 **Ready to Get Started?**
    
    Choose your preferred mode from the sidebar:
    """)
    
    col1, col2, col3 = st.columns([1, 2, 1])
    
    with col2:
        if st.button("🔬 Start with Generator Mode", type="primary", use_container_width=True):
            st.session_state.app_mode = 'Generator'
            st.rerun()
        
        if st.button("📁 Start with Upload Mode", use_container_width=True):
            st.session_state.app_mode = 'Upload'
            st.rerun()
    
    # Example workflow
    st.markdown("""
    ### 📋 **Example Workflow**
    
    **Here's a sample workflow using RNMG:**
    
    1. **Choose Generator Mode** for new reaction mechanisms
    2. **Enter reactions** like: `CO(g) + * ⇌ CO*`
    3. **Add each reaction** one by one using the ⇌ arrow
    4. **Generate matrices** automatically
    5. **Create parameter file** with default values
    6. **Perform mass balance** check (A × ν = 0)
    7. **Run equilibrium analysis** for thermodynamics
    8. **Download results** for further modeling
    
    **Or use Upload Mode** to analyze existing matrix files!
    """)

# Generator Mode
elif selected_mode == "Generator":
    # Reaction input section
    st.subheader("📝 Add Chemical Reactions")

    # Show the equilibrium arrow for easy copying
    st.info("**Use this arrow in your reactions:** ⇌")

    # Input field for reactions
    reaction_input = st.text_input(
        "Enter reaction (e.g., CO(g) + * ⇌ CO*):",
        placeholder="CO(g) + * ⇌ CO*",
        help="Use ⇌ for equilibrium reactions. Supported surface sites: *, _, #, @, &"
    )

    # Buttons for adding reactions
    col1, col2, col3 = st.columns([1, 1, 2])

    with col1:
        if st.button("Add Reaction", type="primary"):
            if reaction_input.strip():
                reactants, products, reaction_str = parse_reaction(reaction_input)
                if reactants is not None and products is not None:
                    st.session_state.reactions.append((len(st.session_state.reactions) + 1, reactants, products, reaction_str))
                    update_matrices()
                    st.success(f"✅ Added: {reaction_str}")
                    st.rerun()
            else:
                st.error("Please enter a reaction")

    with col2:
        if st.button("Clear All"):
            st.session_state.reactions = []
            st.session_state.species_list = []
            st.session_state.atomic_matrix = pd.DataFrame()
            st.session_state.stoich_matrix = pd.DataFrame()
            st.session_state.param_data = []
            st.success("✅ All reactions cleared")
            st.rerun()

    # Display current reactions
    if st.session_state.reactions:
        st.subheader("📋 Current Reactions")
        for i, (_, reactants, products, reaction_str) in enumerate(st.session_state.reactions):
            col1, col2 = st.columns([4, 1])
            with col1:
                st.write(f"**r{i+1}:** {reaction_str}")
            with col2:
                if st.button("Remove", key=f"remove_{i}"):
                    st.session_state.reactions.pop(i)
                    update_matrices()
                    st.rerun()

    # Matrix display section
    if not st.session_state.atomic_matrix.empty:
        st.subheader("📊 Generated Matrices")

        # Matrix consistency check
        if not st.session_state.matrix_consistency['consistent']:
            st.error("⚠️ Matrix inconsistency detected! Species ordering differs between matrices.")

        # Display matrices side by side
        col1, col2 = st.columns(2)

        with col1:
            st.write("**Atomic Matrix (A)**")
            st.caption("Shows atomic composition of each species")
            st.dataframe(st.session_state.atomic_matrix, use_container_width=True)

        with col2:
            st.write("**Stoichiometric Matrix (ν)**")
            st.caption("Shows net stoichiometric coefficients for each reaction")
            st.dataframe(st.session_state.stoich_matrix, use_container_width=True)

        # Generate parameters
        st.subheader("⚙️ Parameters")
        st.caption("Generate and edit kinetic parameters for your reaction system")

        if st.button("Generate Default Parameters"):
            st.session_state.param_data = create_default_parameters()
            st.success("✅ Default parameters generated")

        if st.session_state.param_data:
            param_df = pd.DataFrame(st.session_state.param_data)
            edited_params = st.data_editor(
                param_df,
                use_container_width=True,
                num_rows="dynamic",
                key="param_editor"
            )
            st.session_state.param_data = edited_params.to_dict('records')

        # Analysis section for Generator mode
        if not st.session_state.atomic_matrix.empty and not st.session_state.stoich_matrix.empty:
            st.markdown("---")
            params_df = pd.DataFrame(st.session_state.param_data) if st.session_state.param_data else pd.DataFrame()
            perform_analysis(st.session_state.atomic_matrix, st.session_state.stoich_matrix, params_df)

        # Download section
        st.subheader("💾 Download Matrices")
        col1, col2, col3 = st.columns(3)

        with col1:
            if not st.session_state.atomic_matrix.empty:
                csv_atomic = st.session_state.atomic_matrix.to_csv(index=False)
                st.download_button(
                    label="📥 Download Atomic Matrix",
                    data=csv_atomic,
                    file_name="atomic_matrix.csv",
                    mime="text/csv"
                )

        with col2:
            if not st.session_state.stoich_matrix.empty:
                csv_stoich = st.session_state.stoich_matrix.to_csv(index=False)
                st.download_button(
                    label="📥 Download Stoichiometric Matrix",
                    data=csv_stoich,
                    file_name="stoichiometric_matrix.csv",
                    mime="text/csv"
                )

        with col3:
            if st.session_state.param_data:
                param_df = pd.DataFrame(st.session_state.param_data)
                csv_params = param_df.to_csv(index=False)
                st.download_button(
                    label="📥 Download Parameters",
                    data=csv_params,
                    file_name="parameters.csv",
                    mime="text/csv"
                )

# Upload Mode
elif selected_mode == "Upload":
    # File upload section
    st.subheader("📁 Upload Files")
    st.caption("Upload your matrices and parameters to perform various analyses")

    col1, col2, col3 = st.columns(3)

    with col1:
        uploaded_atomic = st.file_uploader("Upload Atomic Matrix", type=['csv'], key="atomic_upload")
        if uploaded_atomic:
            st.session_state.uploaded_files['atomic'] = pd.read_csv(uploaded_atomic)
            st.success("✅ Atomic matrix uploaded")

    with col2:
        uploaded_stoich = st.file_uploader("Upload Stoichiometric Matrix", type=['csv'], key="stoich_upload")
        if uploaded_stoich:
            st.session_state.uploaded_files['stoich'] = pd.read_csv(uploaded_stoich)
            st.success("✅ Stoichiometric matrix uploaded")

    with col3:
        uploaded_params = st.file_uploader("Upload Parameters", type=['csv'], key="params_upload")
        if uploaded_params:
            st.session_state.uploaded_files['params'] = pd.read_csv(uploaded_params)
            st.success("✅ Parameters uploaded")

    # Display uploaded matrices
    if st.session_state.uploaded_files['atomic'] is not None or st.session_state.uploaded_files['stoich'] is not None:
        st.subheader("📊 Uploaded Matrices")

        col1, col2 = st.columns(2)

        with col1:
            if st.session_state.uploaded_files['atomic'] is not None:
                st.write("**Atomic Matrix (A)**")
                st.caption("Shows atomic composition of each species")
                st.dataframe(st.session_state.uploaded_files['atomic'], use_container_width=True)

        with col2:
            if st.session_state.uploaded_files['stoich'] is not None:
                st.write("**Stoichiometric Matrix (ν)**")
                st.caption("Shows net stoichiometric coefficients for each reaction")
                st.dataframe(st.session_state.uploaded_files['stoich'], use_container_width=True)

    # Display parameters if uploaded
    if st.session_state.uploaded_files['params'] is not None:
        st.subheader("⚙️ Uploaded Parameters")
        st.dataframe(st.session_state.uploaded_files['params'], use_container_width=True)

    # Analysis section for Upload mode
    atomic_matrix = st.session_state.uploaded_files['atomic'] if st.session_state.uploaded_files['atomic'] is not None else pd.DataFrame()
    stoich_matrix = st.session_state.uploaded_files['stoich'] if st.session_state.uploaded_files['stoich'] is not None else pd.DataFrame()
    params_df = st.session_state.uploaded_files['params'] if st.session_state.uploaded_files['params'] is not None else pd.DataFrame()

    if not atomic_matrix.empty or not stoich_matrix.empty:
        st.markdown("---")
        perform_analysis(atomic_matrix, stoich_matrix, params_df)

# Footer
st.markdown("---")
st.markdown("**RNMG v1.1** - Reaction Network Matrix Generator | Built with Streamlit by Kenneth Kusima (tinyurl.com/kennethkusima)")

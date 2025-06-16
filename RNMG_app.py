import streamlit as st
import pandas as pd
import numpy as np
import re
from io import StringIO
from collections import defaultdict, Counter

# Page configuration
st.set_page_config(
    page_title="Reaction Network Matrix Generator",
    page_icon="🔗",
    layout="wide"
)

st.title("Reaction Network Matrix Generator")
st.markdown("Create stoichiometric and atomic matrices from chemical reaction mechanisms")

# Initialize session state
if 'reactions' not in st.session_state:
    st.session_state.reactions = []
if 'species_list' not in st.session_state:
    st.session_state.species_list = []
if 'atomic_matrix' not in st.session_state:
    st.session_state.atomic_matrix = pd.DataFrame()
if 'stoich_matrix' not in st.session_state:
    st.session_state.stoich_matrix = pd.DataFrame()
if 'param_data' not in st.session_state:
    st.session_state.param_data = []
if 'matrix_consistency' not in st.session_state:
    st.session_state.matrix_consistency = {'consistent': True, 'atomic_species': [], 'stoich_species': []}
if 'mass_balance_result' not in st.session_state:
    st.session_state.mass_balance_result = {'balanced': None, 'message': ''}

def parse_species(species_str):
    """Parse a species string to extract species name and stoichiometric coefficient"""
    species_str = species_str.strip()
    
    # Handle coefficients like 2O*, 2*
    coeff_match = re.match(r'^(\d+)(.+)$', species_str)
    if coeff_match:
        coeff = int(coeff_match.group(1))
        species = coeff_match.group(2)
    else:
        coeff = 1
        species = species_str
    
    return species, coeff

def parse_reaction(reaction_str):
    """Parse a reaction string and extract reactants and products"""
    # Replace different arrow types with standard arrow
    reaction_str = re.sub(r'[⇌↔⇄]', '⇌', reaction_str)
    reaction_str = re.sub(r'->', '⇌', reaction_str)
    reaction_str = re.sub(r'<->', '⇌', reaction_str)
    
    if '⇌' not in reaction_str:
        st.error(f"Invalid reaction format: {reaction_str}")
        return None, None, None
    
    # Split into reactants and products
    parts = reaction_str.split('⇌')
    if len(parts) != 2:
        st.error(f"Invalid reaction format: {reaction_str}")
        return None, None, None
    
    reactants_str, products_str = parts
    
    # Parse reactants
    reactants = {}
    for reactant in reactants_str.split('+'):
        reactant = reactant.strip()
        if reactant:
            species, coeff = parse_species(reactant)
            reactants[species] = reactants.get(species, 0) + coeff
    
    # Parse products
    products = {}
    for product in products_str.split('+'):
        product = product.strip()
        if product:
            species, coeff = parse_species(product)
            products[species] = products.get(species, 0) + coeff
    
    return reactants, products, reaction_str

def get_atomic_composition():
    """Define atomic composition for common species"""
    compositions = {
        'CO': {'C': 1, 'O': 1},
        'O2': {'O': 2},
        'CO2': {'C': 1, 'O': 2},
        'CO*': {'C': 1, 'O': 1, '*': 1},
        'O*': {'O': 1, '*': 1},
        'O2*': {'O': 2, '*': 1},
        '*': {'*': 1},
        'H2': {'H': 2},
        'H2O': {'H': 2, 'O': 1},
        'CH4': {'C': 1, 'H': 4},
        'NH3': {'N': 1, 'H': 3},
        'N2': {'N': 2},
        'H*': {'H': 1, '*': 1},
        'OH*': {'O': 1, 'H': 1, '*': 1},
        'C*': {'C': 1, '*': 1},
        'N*': {'N': 1, '*': 1}
    }
    return compositions

def get_species_atomic_composition(species):
    """Get atomic composition for a species, handling gas phase notation"""
    compositions = get_atomic_composition()
    
    # Remove gas phase notation (g), (l), (s) if present
    clean_species = re.sub(r'\([gls]\)', '', species)
    
    # Look up composition
    if clean_species in compositions:
        return compositions[clean_species]
    else:
        # If not found, try to parse simple species
        # This is a basic parser - can be expanded for more complex species
        return parse_species_formula(clean_species)

def parse_species_formula(formula):
    """Basic parser for simple chemical formulas"""
    composition = {}
    
    # Handle surface sites
    if formula == '*':
        return {'*': 1}
    
    # Handle adsorbed species (ending with *)
    if formula.endswith('*'):
        base_formula = formula[:-1]
        composition['*'] = 1
        # Parse the base formula
        base_comp = parse_species_formula(base_formula)
        composition.update(base_comp)
        return composition
    
    # Simple regex pattern to match element-number pairs
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)
    
    for element, count in matches:
        count = int(count) if count else 1
        composition[element] = composition.get(element, 0) + count
    
    return composition if composition else {'Unknown': 1}

def update_matrices():
    """Update atomic and stoichiometric matrices based on current reactions"""
    if not st.session_state.reactions:
        st.session_state.atomic_matrix = pd.DataFrame()
        st.session_state.stoich_matrix = pd.DataFrame()
        st.session_state.species_list = []
        return
    
    # Collect all species
    all_species = set()
    for _, reactants, products, _ in st.session_state.reactions:
        all_species.update(reactants.keys())
        all_species.update(products.keys())
    
    # Separate and sort gas and surface species, with empty sites (*) last
    gas_species = sorted([s for s in all_species if not s.endswith('*')])
    surface_species = [s for s in all_species if s.endswith('*')]
    
    # Sort surface species with empty sites (*) at the end
    other_surface = sorted([s for s in surface_species if s != '*'])
    empty_sites = [s for s in surface_species if s == '*']
    surface_species = other_surface + empty_sites
    
    # Create consistent species order: gas species first, then surface species, with * last
    ordered_species = gas_species + surface_species
    st.session_state.species_list = ordered_species
    
    # Create atomic matrix with consistent ordering
    all_atoms = set()
    
    # Collect all atom types using the new function
    for species in ordered_species:
        composition = get_species_atomic_composition(species)
        all_atoms.update(composition.keys())
    
    # Sort atoms with "*" (surface sites) at the end
    other_atoms = sorted([atom for atom in all_atoms if atom != '*'])
    surface_atom = [atom for atom in all_atoms if atom == '*']
    all_atoms = other_atoms + surface_atom
    
    # Build atomic matrix with consistent species ordering
    atomic_data = []
    for atom in all_atoms:
        row = [atom]
        for species in ordered_species:  # Use ordered_species to ensure consistent ordering
            composition = get_species_atomic_composition(species)
            row.append(composition.get(atom, 0))
        atomic_data.append(row)
    
    columns = ['A\\S'] + ordered_species
    st.session_state.atomic_matrix = pd.DataFrame(atomic_data, columns=columns)
    
    # Create stoichiometric matrix with same species ordering
    stoich_data = []
    
    # Get gas and surface species in the same order as atomic matrix
    gas_species = [s for s in ordered_species if not s.endswith('*')]
    surface_species = [s for s in ordered_species if s.endswith('*')]
    
    # Column headers: P_ for gas phase, theta_ for surface species (with * last)
    # Use the same ordering as atomic matrix
    stoich_columns = ['r\\S'] + [f'P_{s}' for s in gas_species] + [f'theta_{s}' for s in surface_species]
    
    for i, (_, reactants, products, _) in enumerate(st.session_state.reactions):
        row = [f'r{i+1}']
        
        # Gas phase species (in same order as atomic matrix)
        for species in gas_species:
            net_coeff = products.get(species, 0) - reactants.get(species, 0)
            row.append(net_coeff)
        
        # Surface species (in same order as atomic matrix)
        for species in surface_species:
            net_coeff = products.get(species, 0) - reactants.get(species, 0)
            row.append(net_coeff)
        
        stoich_data.append(row)
    
    st.session_state.stoich_matrix = pd.DataFrame(stoich_data, columns=stoich_columns)
    
    # Verification check: ensure both matrices have consistent species ordering
    verify_matrix_consistency()

def verify_matrix_consistency():
    """Verify that atomic and stoichiometric matrices have consistent species ordering"""
    if st.session_state.atomic_matrix.empty or st.session_state.stoich_matrix.empty:
        return True
    
    # Get species from atomic matrix (excluding first column A\S)
    atomic_species = st.session_state.atomic_matrix.columns[1:].tolist()
    
    # Get species from stoichiometric matrix (excluding first column r\S)
    # Remove prefixes P_ and theta_ to get actual species names
    stoich_species = []
    for col in st.session_state.stoich_matrix.columns[1:]:
        if col.startswith('P_'):
            stoich_species.append(col[2:])  # Remove 'P_' prefix
        elif col.startswith('theta_'):
            stoich_species.append(col[6:])  # Remove 'theta_' prefix
    
    # Check if species lists match
    consistency_check = atomic_species == stoich_species
    
    # Store consistency status in session state for display
    st.session_state.matrix_consistency = {
        'consistent': consistency_check,
        'atomic_species': atomic_species,
        'stoich_species': stoich_species
    }
    
    return consistency_check

def check_mass_balance():
    """Check if mass is conserved across all reactions using A*v=0"""
    if st.session_state.atomic_matrix.empty or st.session_state.stoich_matrix.empty:
        return False, "Matrices not available for mass balance check"
    
    # Extract atomic matrix (excluding first column A\S)
    atomic_matrix = st.session_state.atomic_matrix.iloc[:, 1:].values
    
    # Extract stoichiometric matrix (excluding first column r\S)  
    stoich_matrix = st.session_state.stoich_matrix.iloc[:, 1:].values
    
    errors = []
    mass_balanced = True
    
    # Check mass balance for each reaction
    for i, reaction_row in enumerate(stoich_matrix):
        # Compute A * v for reaction i
        mass_balance_result = np.dot(atomic_matrix, reaction_row)
        
        # Check if all elements are approximately zero (mass balanced)
        if not np.allclose(mass_balance_result, 0, atol=1e-10):
            mass_balanced = False
            reaction_name = st.session_state.stoich_matrix.iloc[i, 0]
            errors.append({
                'reaction': i + 1,
                'reaction_name': reaction_name,
                'imbalance': mass_balance_result.tolist()
            })
    
    if mass_balanced:
        return True, "✅ Mass is conserved across all reactions."
    else:
        error_message = "❌ Mass is NOT conserved in the following reactions:\n"
        for error in errors:
            error_message += f"• Reaction {error['reaction']} ({error['reaction_name']}): Imbalance = {error['imbalance']}\n"
        error_message += "\nPlease check and correct the Atomic and/or Stoichiometric matrices."
        return False, error_message

def create_default_parameters():
    """Create default parameter entries for the current reaction set"""
    params = []
    
    # Global parameters
    params.append({'Reaction_Descrp': '', 'Parameter': 'T', 'Values': 320.0, 'Units': 'K'})
    params.append({'Reaction_Descrp': '', 'Parameter': 'R', 'Values': 8.31446, 'Units': 'JK^-1mol^-1'})
    
    # Pressure parameters for gas species (in consistent order)
    gas_species = [s for s in st.session_state.species_list if not s.endswith('*')]
    for i, species in enumerate(gas_species):
        params.append({'Reaction_Descrp': species, 'Parameter': f'P{i+1}', 'Values': 1.0e-8, 'Units': 'bar'})
    
    # Rate constants for reactions
    for i in range(len(st.session_state.reactions)):
        reaction_id = f'r{i+1}'
        params.append({'Reaction_Descrp': reaction_id, 'Parameter': f'k{i+1}f', 'Values': 1.0, 'Units': '-'})
        params.append({'Reaction_Descrp': '', 'Parameter': f'k{i+1}r', 'Values': 1.0, 'Units': '-'})
    
    # Constants
    for i in range(len(st.session_state.reactions)):
        for direction in ['f', 'r']:
            for const_type in ['a', 'b', 'c']:
                params.append({'Reaction_Descrp': 'const', 'Parameter': f'{const_type}{i+1}{direction}', 'Values': 1.0, 'Units': '-'})
    
    return params

# Sidebar for reaction input
with st.sidebar:
    st.header("Reaction Input")
    
    reaction_input = st.text_area(
        "Enter reaction (e.g., CO(g) + * ⇌ CO*):",
        placeholder="CO(g) + * ⇌ CO*",
        help="Use ⇌ for arrows. Use * for surface sites."
    )
    
    if st.button("Add Reaction", type="primary"):
        if reaction_input.strip():
            reactants, products, original = parse_reaction(reaction_input.strip())
            if reactants is not None and products is not None:
                st.session_state.reactions.append((len(st.session_state.reactions), reactants, products, original))
                update_matrices()
                st.session_state.param_data = create_default_parameters()
                st.session_state.mass_balance_result = {'balanced': None, 'message': ''}  # Reset mass balance
                st.success("Reaction added!")
                st.rerun()
    
    if st.session_state.reactions:
        st.subheader("Current Reactions")
        for i, (idx, reactants, products, original) in enumerate(st.session_state.reactions):
            col1, col2 = st.columns([3, 1])
            with col1:
                st.write(f"r{i+1}: {original}")
            with col2:
                if st.button("❌", key=f"del_{i}", help="Delete reaction"):
                    st.session_state.reactions.pop(i)
                    update_matrices()
                    st.session_state.param_data = create_default_parameters()
                    st.session_state.mass_balance_result = {'balanced': None, 'message': ''}  # Reset mass balance
                    st.rerun()
    
    if st.button("Clear All", type="secondary"):
        st.session_state.reactions = []
        st.session_state.species_list = []
        st.session_state.atomic_matrix = pd.DataFrame()
        st.session_state.stoich_matrix = pd.DataFrame()
        st.session_state.param_data = []
        st.session_state.mass_balance_result = {'balanced': None, 'message': ''}
        st.rerun()

# Main content area
if st.session_state.reactions:
    # Tabs for different matrices
    tab1, tab2, tab3 = st.tabs(["📊 Stoichiometric Matrix", "⚛️ Atomic Matrix", "⚙️ Parameters"])
    
    with tab1:
        st.subheader("Stoichiometric Matrix")
        if not st.session_state.stoich_matrix.empty:
            st.dataframe(st.session_state.stoich_matrix, use_container_width=True)
            
            # Download section
            col1, col2 = st.columns([2, 1])
            with col1:
                stoich_filename = st.text_input("Stoichiometric matrix filename:", value="Stoich_1.csv")
            with col2:
                csv_stoich = st.session_state.stoich_matrix.to_csv(index=False)
                st.download_button(
                    label="📥 Download CSV",
                    data=csv_stoich,
                    file_name=stoich_filename,
                    mime="text/csv"
                )
    
    with tab2:
        st.subheader("Atomic Matrix")
        if not st.session_state.atomic_matrix.empty:
            st.dataframe(st.session_state.atomic_matrix, use_container_width=True)
            
            # Download section
            col1, col2 = st.columns([2, 1])
            with col1:
                atomic_filename = st.text_input("Atomic matrix filename:", value="Atomic_1.csv")
            with col2:
                csv_atomic = st.session_state.atomic_matrix.to_csv(index=False)
                st.download_button(
                    label="📥 Download CSV",
                    data=csv_atomic,
                    file_name=atomic_filename,
                    mime="text/csv"
                )
    
    with tab3:
        st.subheader("Parameter Table")
        
        if st.session_state.param_data:
            # Create editable dataframe
            param_df = pd.DataFrame(st.session_state.param_data)
            
            # Display editable table
            edited_df = st.data_editor(
                param_df,
                use_container_width=True,
                num_rows="dynamic",
                column_config={
                    "Reaction_Descrp": st.column_config.TextColumn("Reaction Description"),
                    "Parameter": st.column_config.TextColumn("Parameter"),
                    "Values": st.column_config.NumberColumn("Values", format="%.2e"),
                    "Units": st.column_config.TextColumn("Units")
                }
            )
            
            # Update session state with edits
            st.session_state.param_data = edited_df.to_dict('records')
            
            # Download section
            col1, col2 = st.columns([2, 1])
            with col1:
                param_filename = st.text_input("Parameter file filename:", value="Param_1.csv")
            with col2:
                csv_param = edited_df.to_csv(index=False)
                st.download_button(
                    label="📥 Download CSV",
                    data=csv_param,
                    file_name=param_filename,
                    mime="text/csv"
                )
    
    # Matrix consistency check
    st.markdown("---")
    if 'matrix_consistency' in st.session_state:
        if st.session_state.matrix_consistency['consistent']:
            st.success("✅ Matrix Consistency Check: Both atomic and stoichiometric matrices have consistent species ordering (Gas → Surface)")
        else:
            st.error("❌ Matrix Consistency Check: Species ordering mismatch detected!")
            with st.expander("View Ordering Details"):
                col1, col2 = st.columns(2)
                with col1:
                    st.write("**Atomic Matrix Species Order:**")
                    for i, species in enumerate(st.session_state.matrix_consistency['atomic_species']):
                        st.write(f"{i+1}. {species}")
                with col2:
                    st.write("**Stoichiometric Matrix Species Order:**")
                    for i, species in enumerate(st.session_state.matrix_consistency['stoich_species']):
                        st.write(f"{i+1}. {species}")
    
    # Mass Balance Check Section
    st.markdown("### 🧮 Mass Balance Verification")
    col1, col2 = st.columns([1, 3])
    
    with col1:
        if st.button("Check Mass Balance", type="primary"):
            balanced, message = check_mass_balance()
            st.session_state.mass_balance_result = {'balanced': balanced, 'message': message}
            st.rerun()
    
    with col2:
        if st.session_state.mass_balance_result['balanced'] is not None:
            if st.session_state.mass_balance_result['balanced']:
                st.success(st.session_state.mass_balance_result['message'])
            else:
                st.error(st.session_state.mass_balance_result['message'])
    
    # Information about mass balance
    with st.expander("ℹ️ About Mass Balance Check"):
        st.markdown("""
        **Mass Balance Verification** ensures that atoms are conserved across all reactions.
        
        **How it works:**
        - For each reaction, the atomic matrix (A) is multiplied by the stoichiometric coefficients (v)
        - If A × v = 0 for all reactions, mass is conserved
        - Any non-zero result indicates a mass imbalance
        
        **Example:**
        For reaction: CO + O* → CO2 + *
        - Carbon: 1×1 + 0×1 = 1 (reactants) vs 1×1 + 0×1 = 1 (products) ✓
        - Oxygen: 1×1 + 1×1 = 2 (reactants) vs 2×1 + 0×1 = 2 (products) ✓
        """)
    
    # Summary section
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Reactions", len(st.session_state.reactions))
    with col2:
        st.metric("Species", len(st.session_state.species_list))
    with col3:
        gas_count = len([s for s in st.session_state.species_list if not s.endswith('*')])
        st.metric("Gas Species", gas_count)
    with col4:
        surface_count = len([s for s in st.session_state.species_list if s.endswith('*')])
        st.metric("Surface Species", surface_count)
    
    # Species ordering display
    with st.expander("📋 Current Species Ordering (Gas → Surface, * Last)"):
        col1, col2 = st.columns(2)
        gas_species = [s for s in st.session_state.species_list if not s.endswith('*')]
        surface_species = [s for s in st.session_state.species_list if s.endswith('*')]
        
        with col1:
            st.write("**Gas Species:**")
            for i, species in enumerate(gas_species):
                st.write(f"{i+1}. {species}")
        
        with col2:
            st.write("**Surface Species (Empty Sites Last):**")
            for i, species in enumerate(surface_species):
                if species == '*':
                    st.write(f"{i+1}. {species} ← Empty Sites")
                else:
                    st.write(f"{i+1}. {species}")

else:
    # Welcome message and example
    st.markdown("""
    ## 🚀 Getting Started
    
    Enter chemical reactions in the sidebar to automatically generate:
    - **Stoichiometric Matrix**: Shows how each reaction affects each species
    - **Atomic Matrix**: Shows atomic composition of each species  
    - **Parameter Table**: Editable parameters for reaction kinetics
    
    All matrices maintain consistent ordering: **Gas → Surface → Empty Sites (*)** last.
    
    ### 📝 Example Reactions:
    ```
    CO(g) + * ⇌ CO*
    O2(g) + * ⇌ O2*
    O2* + * ⇌ 2O*
    CO* + O* ⇌ CO2(g) + 2*
    ```
    
    ### 💡 Tips:
    - Use `*` to represent surface sites
    - Use `(g)` to specify gas phase (optional)
    - Supported arrows: `⇌`
    - Add coefficients like `2O*` for multiple species
    - **Automatic Ordering**: Gas species columns appear first, then surface species
    - **Empty Sites Last**: `*` (empty sites) always appear in the last column/row
    - **Consistency Check**: Both matrices maintain identical species ordering
    - **Mass Balance**: Use the mass balance checker to verify atom conservation
    """)

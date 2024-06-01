# Define a function to print in a tabular format
def print_configs(config_list):
    # Get all keys for header
    headers = list(config_list[0].keys())
    
    # Calculate the width of each column
    col_widths = {key: max(len(str(key)), max(len(str(config[key])) for config in config_list)) for key in headers}
    
    # Print headers
    header_line = ' | '.join(f"{key:<{col_widths[key]}}" for key in headers)
    separator = '-+-'.join('-' * col_widths[key] for key in headers)
    print(header_line)
    print(separator)
    
    # Print each configuration
    for config in config_list:
        line = ' | '.join(f"{str(config[key]):<{col_widths[key]}}" for key in headers)
        print(line)

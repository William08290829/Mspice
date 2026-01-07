import sys
import os
import uuid
import io
import base64
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from flask import Flask, request, jsonify, render_template, send_from_directory

# --- CONFIG ---
BACKEND_DIR = os.path.join(os.path.dirname(__file__), 'Python_Based_SPICE_Circuit_Simulation', 'src')
sys.path.append(BACKEND_DIR)

# --- MOCKS ---
# Mock utils to capture plot instead of showing it
import utils

last_plot_b64 = None
last_csv_data = None

def mock_plot_picture(templist, namelist, times, result, cmd, x, circuit_name):
    global last_plot_b64
    print("Mock Plot Called")
    
    plt.figure(figsize=(10, 6))
    
    # Reconstruct node mapping
    dic_result={}
    nodes=[]
    counter=1
    for i in range(len(templist)):
        if(templist[i]=="0"):
            dic_result["gnd"]=0
            # nodes.append("gnd")
        else:
            dic_result["v("+templist[i]+")"]=counter
            nodes.append("v("+templist[i]+")")
            counter=counter+1
            
    # Plot all voltage nodes (simple auto-plot)
    for node in nodes:
        idx = dic_result[node] - 1
        y = []
        for r in result:
             # Depending on AC/DC/Tran, r might be complex or float
             val = r[idx]
             if isinstance(val, complex) or isinstance(val, np.complex128):
                y.append(abs(val)) # Magnitude for AC
             else:
                y.append(val)
        
        plt.plot(x, y, label=node)

    plt.title(circuit_name)
    plt.grid(True)
    plt.legend()
    
    # AC scale handling
    if cmd[0] == "ac" and (len(cmd)>1 and (cmd[1] == "dec" or cmd[1] == "oct")):
        plt.xscale("log")
        plt.xlabel("Frequency (Hz)")
    elif cmd[0] == "dc":
        plt.xlabel("Sweep Voltage (V)")
    else:
        plt.xlabel("Time/Frequency")

    plt.ylabel("Magnitude / Voltage")

    # Save to buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    last_plot_b64 = base64.b64encode(buf.getvalue()).decode('utf-8')
    plt.close()

    # --- Capture Data for CSV ---
    # Construct headers
    # Determine X label based on command
    x_label = "Time (s)"
    if cmd[0] == "ac": x_label = "Frequency (Hz)"
    elif cmd[0] == "dc": x_label = "Sweep Voltage (V)"
    
    # Add units to headers (Assuming Voltage for now)
    headers = [x_label] + [f"{n} (V)" for n in nodes]
    
    # Construct rows
    # Gather Y-values
    # Note: re-looping to maintain consistent order with headers
    y_lists = [] # List of lists [[y1_points...], [y2_points...]]
    for node in nodes:
        idx = dic_result[node] - 1
        y_val_list = []
        for r in result:
             val = r[idx]
             if isinstance(val, complex) or isinstance(val, np.complex128):
                y_val_list.append(abs(val))
             else:
                y_val_list.append(val)
        y_lists.append(y_val_list)
    
    # Transpose to rows: [x1, y1_node1, y1_node2...], [x2, y2_node1, y2_node2...]
    rows = []
    for i in range(len(x)):
        row = [x[i]]
        for y_l in y_lists:
            row.append(y_l[i])
        rows.append(row)

    global last_csv_data
    last_csv_data = {"headers": headers, "rows": rows}

# Override
utils.plot_picture = mock_plot_picture
utils.plot_picture_trans = mock_plot_picture
utils.plot_picture_dc = mock_plot_picture # Assuming signature matches or similar enough for simple plot
# Note: plot_picture_dc has extra arg 'object_source', we might need a separate mock if signature differs too much
# Checking `plot_picture_dc` signature in utils.py: 
# def plot_picture_dc(templist,namelist,times,result,cmd,x,circuit_name,object_source):
# Our mock above doesn't have object_source. We need to handle it.

def mock_plot_picture_dc(templist, namelist, times, result, cmd, x, circuit_name, object_source):
    # Just call the main mock ignoring object_source for now, or pass it
    mock_plot_picture(templist, namelist, times, result, cmd, x, circuit_name)

def mock_plot_picture_trans(templist, namelist, result, cmd, x, circuit_name):
    # Transient analysis doesn't pass 'times' explicitly in the same slot, but 'x' is time points
    # We can infer 'times' as len(x) or len(result)
    mock_plot_picture(templist, namelist, len(result), result, cmd, x, circuit_name)

utils.plot_picture_dc = mock_plot_picture_dc
utils.plot_picture_trans = mock_plot_picture_trans


# --- IMPORTS AFTER MOCK ---
from new_parser import parse
from calculation_module_0129 import ac_analysis, Dc_analysis
from dc_bias import Dc_bias
from transient_0327 import transient_analysis

# Configure Flask to serve static files from 'static' folder
# static_url_path='/static' means requests to /static/... go to static_folder
app = Flask(__name__, static_url_path='/static', static_folder='static', template_folder='.')

@app.route('/')
def home():
    # Landing Page
    return render_template('index.html')

@app.route('/index.html')
def home_explicit():
    return render_template('index.html')

@app.route('/index-zh.html')
def home_zh():
    return render_template('index-zh.html')

@app.route('/app')
def app_page():
    # App Interface
    return render_template('app.html')

@app.route('/api/simulate', methods=['POST'])
def simulate():
    global last_plot_b64, last_csv_data
    last_plot_b64 = None
    last_csv_data = None
    
    data = request.json
    netlist_text = data.get('netlist')
    if not netlist_text:
        return jsonify({"error": "No netlist provided"}), 400

    cmd = None
    cleaned_netlist_lines = []
    
    # Parse lines to find command and strip it from file
    lines = netlist_text.splitlines()
    for line in lines:
        stripped = line.strip().lower()
        if stripped.startswith(".ac") or stripped.startswith(".dc") or stripped.startswith(".tran"):
            cmd = stripped[1:].split() # Found command
            # Do NOT add to cleaned_netlist_lines (or comment it out) to prevent parser error
            cleaned_netlist_lines.append("* " + line) # Comment out command in file
        else:
            cleaned_netlist_lines.append(line)

    if not cmd:
        return jsonify({"error": "No simulation command (.ac, .dc) found in netlist."}), 400

    # Generate unique filename for this request
    unique_filename = f"temp_{uuid.uuid4().hex}.sp"

    try:
        # Write temp file for parser (without active command lines)
        with open(unique_filename, "w") as f:
            f.write("\n".join(cleaned_netlist_lines))

        # --- LOGIC PORTED FROM test_func.py ---
        circuit_name, components, namelist = parse(unique_filename)
        if circuit_name.startswith('*'):
             circuit_name = circuit_name.lstrip('*').strip()


        # Setup Matrix Logic (Copied from test_func.py)
        typelist=[c.type_ for c in components]
        nodelist = [c.nodes for c in components]
        valuelist = [c.value for c in components]
        element=len(namelist)
        templist=[]
        for i in range(element):
            for j in range(2):
                templist.append(nodelist[i][j])
        templist=list(dict.fromkeys(templist))
        
        matrix_size=0
        for i in range(len(templist)):
            if(templist[i]=="0"):
                continue
            elif(templist[i]!=""):
                matrix_size=matrix_size+1
            else:
                continue
        element2=matrix_size
        
        for i in range(element):
            if(namelist[i][0]=="v"):
                matrix_size=matrix_size+1
            elif(namelist[i][0]=="l"):
                matrix_size=matrix_size+1

        dic_node={}
        counter=1
        for i in range(len(templist)):
            if(templist[i]=='0'):
                dic_node["0"]=0
            else:
                temp_string=templist[i]
                dic_node[temp_string]=counter
                counter=counter+1

        # EXECUTE COMMAND
        if cmd[0] == "ac":
            # Parsing args
            # cmd template: ac dec 100 1 100k
            # test_func.py logic:
            start_freq = utils.unit_symbol(cmd[-2])
            end_freq = utils.unit_symbol(cmd[-1])
            freq_step = utils.unit_symbol(cmd[-3])
            
            vds_his=[]
            vgs_his=[]
            
            # AC Analysis
            result=[]
            place=[]
            dic_result={}
            
            ac_analysis(
                templist, place, element2, element, dic_node, matrix_size, cmd,
                result, namelist, nodelist, circuit_name, valuelist, typelist, components,
                dic_result, freq_step, start_freq, end_freq, vds_his, vgs_his
            )
            
        elif cmd[0] == "dc":
            # dc Vin 0 5 0.1
            start_V = utils.unit_symbol(cmd[-3])
            end_V = utils.unit_symbol(cmd[-2])
            V_step = utils.unit_symbol(cmd[-1])
            object_source = cmd[1]
            
            vds_his=[]
            vgs_his=[]
            null=[]
            C_pass=[[]]
            L_pass=[[]]
            
            # DC Bias (Initial point?) - test_func.py calls Dc_bias first
            # We skip specific history tracking for simplicity unless strictly needed for non-linear?
            # test_func.py lines 169-189 call Dc_bias
            
            nonlin_his,vds_his,vgs_his,C_pass,L_pass,X0 = Dc_bias(
                templist, element2, element, dic_node, matrix_size, cmd, namelist, nodelist,
                circuit_name, valuelist, typelist, components, object_source,
                vgs_his, vds_his, null, start_V, C_pass, L_pass
            )
            
            Dc_analysis(
               templist, element2, element, dic_node, matrix_size, cmd, namelist, nodelist,
               circuit_name, valuelist, typelist, components, V_step, start_V, end_V,
               object_source, vds_his, vgs_his
            )

        elif cmd[0] == "tran":
            # tran command logic
            # Call the refactored transient_analysis function
            # Argument matches the function signature we defined
            transient_analysis(circuit_name, components, namelist, cmd)

        else:
             return jsonify({"error": f"Command {cmd[0]} not supported yet."}), 400

        if last_plot_b64:
            return jsonify({"image": last_plot_b64, "data": last_csv_data})
        else:
            return jsonify({"message": "Simulation ran but no graph produced."})

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500

    finally:
        # Cleanup: Delete the temp file to prevent disk clutter
        if os.path.exists(unique_filename):
            try:
                os.remove(unique_filename)
            except Exception as cleanup_error:
                print(f"Failed to delete temp file {unique_filename}: {cleanup_error}")

if __name__ == '__main__':
    print("Starting Mspice Server...")
    app.run(port=5000, debug=True)

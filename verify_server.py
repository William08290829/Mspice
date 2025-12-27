
import requests
import base64

url = "http://127.0.0.1:5000/api/simulate"

# Sample Netlist from testcase (or simple one)
# Need a valid netlist to run.
# Payload with Netlist
payload = {
    "netlist": """
* RC Transient Test
V1 1 0 PULSE(0 5 1m 1n 1n 1m 2m)
R1 1 2 1k
C1 2 0 1u
.tran 10u 5m
.end
"""
}

try:
    print("Sending Request to", url)
    response = requests.post(url, json=payload)
    response.raise_for_status()
    
    data = response.json()
    if "image" in data:
        print("Success! Image received.")
        print("Image length:", len(data["image"]))
        # Save to see if it's valid
        with open("test_output.png", "wb") as f:
            f.write(base64.b64decode(data["image"]))
        print("Saved test_output.png")
    else:
        print("Response OK but no image:", data)

except Exception as e:
    print("Verification Failed:", e)
    if 'response' in locals():
        print("Response Text:", response.text)

import numpy as np
import pyproj
import matplotlib.pyplot as plt
def convert_rd_to_wgs84(points_rd):
    # Define the RD New (Amersfoort) CRS
    proj_rd_new = pyproj.Proj(init="epsg:28992")

    # Define the WGS84 CRS (used for lat/lon)
    proj_wgs84 = pyproj.Proj(init="epsg:4326")

    # Convert RD New to WGS84 (Lat, Lon) and retain height
    points_wgs84 = [
        (*pyproj.transform(proj_rd_new, proj_wgs84, point[0], point[1]), point[2])
        for point in points_rd
    ]
    
    return points_wgs84

from PIL import Image

def get_image_coordinates(img_path, center_x, center_y, block_size=10, color=(255, 0, 0)):
    """
    Highlights a square area of specified size around a center pixel in the given color.
    
    :param img_path: Path to the image file
    :param center_x: X coordinate of the center pixel
    :param center_y: Y coordinate of the center pixel
    :param block_size: Size of the square block to color (default is 10x10)
    :param color: RGB color to fill the block (default is red)
    """
    # Open the image
    img = Image.open(img_path)
    
    # Calculate top-left corner of the block
    pix_x = center_x - block_size // 2
    pix_y = center_y - block_size // 2
    
    # Modify a block_size x block_size area around the center pixel
    for i in range(pix_x, pix_x + block_size):
        for j in range(pix_y, pix_y + block_size):
            img.putpixel((i, j), color)
    
    # Show the modified image
    img.show()
    # Optionally save or display the image with matplotlib
    # plt.imshow(img)

# Example usage




def compute_projection_matrix(points_rd_array, image_points):
    """
    Computes the projection matrix P using real-world coordinates and image coordinates.
    
    :param points_rd_array: Array of real-world coordinates, shape (N, 3) where N is the number of points
    :param image_points: Array of image coordinates, shape (N, 2)
    :return: Projection matrix P of shape (3, 4)
    """
    # Number of points
    N = points_rd_array.shape[0]
    
    # Initialize A matrix and b vector
    A = []
    b = []
    
    # Construct the A matrix and b vector based on the provided points
    for i in range(N):
        X, Y, Z = points_rd_array[i, :]
        x_prime, y_prime = image_points[i, :]
        
        # Append rows for the A matrix
        A.append([X, Y, Z, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        A.append([0, 0, 0, 0, X, Y, Z, 1, 0, 0, 0, 0])
        A.append([0, 0, 0, 0, 0, 0, 0, 0, X, Y, Z, 1])
        
        # Append corresponding values for the b vector
        b.append(x_prime)
        b.append(y_prime)
        b.append(1)
    
    # Convert A and b to numpy arrays
    A = np.array(A)
    b = np.array(b)
    
    # Solve the least-squares problem to find the projection parameters p
    p, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    
    # Reshape the solution into a 3x4 matrix
    P = p.reshape(3, 4)
    
    return P

def plot_gcp_residuals(gcp_array, residuals, padding_x, padding_y, title="GCP Residuals Plot", x_label="Easting direction", y_label="Northing direction", scale=2, resolution=(8, 6)):
    """
    Plots the Ground Control Points (GCP) residuals along with the Argus Tower position on a quiver plot.
    
    :param gcp_array: Array of GCP coordinates (easting, northing, elevation), shape (N, 3)
    :param residuals: Array of residual errors (easting_error, northing_error), shape (N, 2)
    :param argus_coordinates: Tuple with the coordinates of the Argus Tower (easting, northing)
    :param padding_x: Padding added to x-axis limits
    :param padding_y: Padding added to y-axis limits
    :param title: Title of the plot
    :param x_label: Label for the x-axis
    :param y_label: Label for the y-axis
    :param scale: Scale for the residuals in the quiver plot (default is 15)
    :param resolution: Resolution of the plot figure in inches (width, height), default (8, 6)
    :return: Displays the plot and outputs the image
    """
    # Extract easting and northing components for GCPs
    easting = gcp_array[:, 0]
    northing = gcp_array[:, 1]

    # Extract residual errors for easting and northing
    easting_errors = residuals[:, 0]
    northing_errors = residuals[:, 1]

    # Define GCP labels
    gcp_labels = ["GCP4", "GCP5", "GCP6", "GCP7", "GCP8", "GCP9", "GCP10"]

    # Calculate axis limits with padding
    min_easting, max_easting = np.min(easting) - padding_x, np.max(easting) + padding_x
    min_northing, max_northing = np.min(northing) - padding_y, np.max(northing) + padding_y

    # Define offset for the labels
    x_offset = 6
    y_offset = 6



    # Create a quiver plot
    plt.figure(figsize=resolution)
    plt.quiver(easting, northing, easting_errors, northing_errors, angles='xy', scale_units='xy', scale=scale, color='blue')
    plt.scatter(easting, northing, color='red', label='GCPs')  # Mark GCP locations

    # Set axis limits with padding
    plt.xlim(min_easting, max_easting)
    plt.ylim(min_northing, max_northing)

    # Flip the y-axis and x-axis if needed
    plt.gca().invert_yaxis()
   

    # Add labels for each GCP with offset
    for i, label in enumerate(gcp_labels):
        plt.text(easting[i] + x_offset, northing[i] + y_offset, label, fontsize=10, ha='right', va='bottom')

    # Add labels and grid
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend()
    plt.grid(True)

    # Display the plot
    plt.show()
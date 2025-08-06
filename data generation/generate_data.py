
import numpy as np
import scipy.stats as stats

def generate_gaussian(n, dim_joint, dim_ind, mu_1, mu_2, sigma_1, sigma_2):
    pdf = [1]*n
    with open(r'/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/test_data.txt', 'w') as fp:
        for i in range(n):
            v_threshold = np.random.uniform()
            if v_threshold < 0.4:
                mean_1, var_1 = mu_1[0], sigma_1[0]
            else:
                mean_1, var_1 = mu_1[1], sigma_1[1]
            x_joint = np.random.multivariate_normal(mean_1, var_1).tolist()
            pdf[i] *= 0.4 * stats.multivariate_normal(mu_1[0], sigma_1[0]).pdf(x_joint) \
                + 0.6 * stats.multivariate_normal(mu_1[1], sigma_1[1]).pdf(x_joint)
            
            x_rest = []
            for j in range(dim_ind):
                if np.random.uniform() < 0.5:
                    mean_2, var_2 = mu_2[0], sigma_2[0]
                else:
                    mean_2, var_2 = mu_2[1], sigma_2[1]
                x_rest.append(np.random.normal(mean_2, var_2))
                pdf[i] *= 0.5*stats.norm(mu_2[0], sigma_2[0]).pdf(x_rest[-1]) \
                    + 0.5*stats.norm(mu_2[1], sigma_2[1]).pdf(x_rest[-1])
            x = x_joint+x_rest
            for value in x:
                # write each item on a new line
                fp.write("%.10f\n" % value)
            print('Done')
            
    with open(r'/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/test_data_pdf.txt', 'w') as fp:
        for value in pdf:
            # write each item on a new line
            fp.write("%.20f\n" % value)
        print('Done 2')

    return


def generate_beta():
    return


if __name__ == "__main__":
    n = 5000
    dim_joint = 5
    dim_ind = 5
    mu_1 = [[0.3]*dim_joint, [0.7]*dim_joint]
    mu_2 = [0.3,0.7]
    sigma_1 = [[[1e-3]*dim_joint for i in range(dim_joint)],[[1e-3]*dim_joint for i in range(dim_joint)]]
    sigma_2 = [0.1,0.1]
    for i in range(dim_joint):
        sigma_1[0][i][i] += 1.5e-3
        sigma_1[1][i][i] += 1.5e-3
    
    generate_gaussian(n, dim_joint, dim_ind, mu_1, mu_2, sigma_1, sigma_2)
    
    


using System;
using System.Runtime.InteropServices;
using System.Security;

namespace pardiso {

	public sealed class Pardiso {
		private Pardiso() {
		}

		public static int pardiso(IntPtr[] handle,
				ref int maxfct, ref int mnum,
				ref int mtype, ref int phase, ref int n,
				double[] a, int[] ia, int[] ja, int[] perm,
				ref int nrhs, int[] iparm, ref int msglvl,
				double[] b, double[] x, ref int error) {
			return PardisoNative.pardiso(handle, ref maxfct, ref mnum, ref mtype,
					ref phase, ref n, a, ia, ja, perm, ref nrhs, iparm, ref msglvl, b, x, ref error);
		}
	}

	[SuppressUnmanagedCodeSecurity]
	internal sealed class PardisoNative {
		private PardisoNative() {
		}

		[DllImport("libpardiso600-WIN-X86-64.dll", CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
		internal static extern int pardiso(
				[In, Out] IntPtr[] handle,
				ref int maxfct, 
				ref int mnum,
				ref int mtype, 
				ref int phase, 
				ref int n,
				[In] double[] a,
				[In] int[] ia,
				[In] int[] ja,
				[In] int[] perm,
				ref int nrhs, 
				[In, Out] int[] iparm, 
				ref int msglvl,
				[In, Out] double[] b, 
				[Out] double[] x, 
				ref int error);
	}
}
